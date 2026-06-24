#include "transport/osh_transport_neutron.h"

#include <math.h>

#include "common/osh_diag.h"
#include "common/osh_step.h"
#include "gemca/runtime/osh_gemca_runtime.h"
#include "material/runtime/osh_material_runtime.h"
#include "openshieldhit/const.h"
#include "openshieldhit/geometry.h"
#include "openshieldhit/material.h"
#include "particle/osh_particle.h"
#include "particle/osh_particle_const.h"
#include "particle/osh_particle_pdg.h"
#include "physics/neutron/osh_neutron_reaction.h"
#include "physics/neutron/osh_neutron_xsec.h"
#include "physics/nuclear/osh_nuclear_handler.h"
#include "random/osh_rng.h"
#include "scoring/runtime/osh_scoring_step.h"
#include "transport/osh_neutron_pool.h"
#include "transport/osh_transport.h"
#include "transport/osh_transport_boundary.h"

/* --------------------------------------------------------------------------
 * Static helpers
 * -------------------------------------------------------------------------- */

static struct particle const k_neutron_species = {OSH_PART_MASS_NEUTRON, OSH_PART_PDG_NEUTRON, 0, 0u, 1u, 0u};

/*
 * Total macroscopic cross section Σ_tot [cm⁻¹] for one material cell.
 * Σ_tot = Σ_i n_i σ_tot,i,  with n_i [cm⁻³] and σ in cm².
 *
 * Called separately from osh_neutron_reaction_sample() so the free-path
 * length can be sampled before committing to a full reaction.  The two
 * xsec_lookup calls per element are redundant but acceptable for the
 * minimal model; they share the same xsec warning-once state.
 */
static double neutron_sigma_tot_cm(struct osh_neutron_xsec *xsec,
                                   struct osh_nuclear_handler const *handler,
                                   size_t mat_idx,
                                   double rho_g_cm3,
                                   double e_mev) {
    struct osh_nuclear_elem const *elems; /* element list for this material */
    size_t n_elems;
    size_t i;
    double sum; /* accumulator in mb·cm⁻³ */
    struct osh_neutron_xsec_result sig;

    elems = handler->elem_pool + handler->elem_offset[mat_idx];
    n_elems = handler->elem_count[mat_idx];
    sum = 0.0;

    for (i = 0u; i < n_elems; ++i) {
        double nd_i; /* number density of element i [cm⁻³] */
        unsigned int a;
        a = osh_neutron_xsec_resolve_a(elems[i].z, elems[i].a);
        nd_i = (double) elems[i].mass_fraction * rho_g_cm3 * OSH_NAVOGADRO / (double) a;
        osh_neutron_xsec_lookup(xsec, (int) elems[i].z, (int) a, e_mev, &sig);
        sum += nd_i * sig.tot;
    }
    return sum * OSH_MB_TO_CM2; /* convert mb·cm⁻³ → cm⁻¹ */
}

/*
 * Advance neutron slot k to the zone boundary plus a small eps nudge so the
 * next zone-ref query resolves the new zone correctly.
 */
static void advance_to_boundary(struct osh_neutron_pool *pool, size_t k, double dist) {
    osh_transport_advance_to_boundary(
        &pool->x[k], &pool->y[k], &pool->z[k], pool->ux[k], pool->uy[k], pool->uz[k], dist);
}

static enum osh_status score_neutron_step(struct osh_scoring_runtime *score_rt,
                                          struct osh_neutron_pool const *pool,
                                          size_t k,
                                          struct osh_zone_ref const *zone_ref,
                                          double rho,
                                          double step_len) {
    struct step st;

    if (!score_rt || !(step_len > 0.0)) {
        return OSH_OK;
    }

    st.p[0] = pool->x[k];
    st.p[1] = pool->y[k];
    st.p[2] = pool->z[k];
    st.p[3] = pool->e[k];
    st.q[0] = pool->x[k] + pool->ux[k] * step_len;
    st.q[1] = pool->y[k] + pool->uy[k] * step_len;
    st.q[2] = pool->z[k] + pool->uz[k] * step_len;
    st.q[3] = pool->e[k];
    st.v[0] = pool->ux[k];
    st.v[1] = pool->uy[k];
    st.v[2] = pool->uz[k];
    st.w[0] = pool->ux[k];
    st.w[1] = pool->uy[k];
    st.w[2] = pool->uz[k];
    st.ds = step_len;
    st.de = 0.0;
    st.rho = rho;
    st.wt = pool->wt[k];
    st.voxel_idx = (zone_ref && zone_ref->has_hu) ? zone_ref->voxel_idx : 0u;
    st.medium = zone_ref ? (int) zone_ref->material_idx : -1;
    st.zone = zone_ref ? (int) zone_ref->zone_idx : -1;
    st.system = OSH_COORD_UNIVERSE;
    st.prim_idx = pool->prim_idx[k];
    st.gen = pool->gen[k];
    st.has_voxel = (zone_ref && zone_ref->has_hu) ? 1u : 0u;

    return osh_scoring_score_step(score_rt, &k_neutron_species, &st);
}

/*
 * Push a neutron secondary from a compound/FBU event back into the pool.
 * Position is inherited from the parent slot; direction and energy come from
 * the secondary descriptor.  RNG stream is split from the parent so the child
 * history is independent from all sibling histories.
 */
static void push_neutron_secondary(struct osh_neutron_pool *pool, size_t k, struct osh_nuclear_secondary const *sec) {
    size_t slot; /* index of the new entry */

    if (pool->n >= pool->capacity) {
        pool->n_dropped++;
        return;
    }
    slot = pool->n++;
    pool->n_created++;
    pool->x[slot] = pool->x[k];
    pool->y[slot] = pool->y[k];
    pool->z[slot] = pool->z[k];
    pool->ux[slot] = sec->dir[0];
    pool->uy[slot] = sec->dir[1];
    pool->uz[slot] = sec->dir[2];
    pool->e[slot] = sec->energy;
    pool->wt[slot] = pool->wt[k];
    pool->prim_idx[slot] = pool->prim_idx[k];
    pool->gen[slot] = (pool->gen[k] < 255u) ? (uint8_t) (pool->gen[k] + 1u) : 255u;
    osh_rng_split(&pool->rng[slot], &pool->rng[k]);
}

/*
 * Apply a reaction event to pool slot k.
 *
 * ELASTIC: update direction and energy in pool; local_deposit_mev (heavy
 *          recoil) is currently discarded — scoring call goes here once
 *          osh_scoring_score_point() is implemented.
 *
 * COMPOUND: neutron secondaries are pushed back to the pool for the next
 *           wavefront pass; ion secondaries would go to the ion pool with
 *           ion-feedback enabled (not yet wired).
 *
 * All other channels: kill the neutron (e = 0); energy goes to
 *   local_deposit_mev which is currently discarded (same scoring TODO).
 */
static void apply_event(struct osh_neutron_pool *pool, size_t k, struct osh_neutron_reaction_event const *ev) {
    size_t i;
    struct osh_nuclear_secondary const *sec;

    switch (ev->kind) {
    case OSH_NEUTRON_REACTION_ELASTIC:
        pool->ux[k] = ev->neutron_dir[0];
        pool->uy[k] = ev->neutron_dir[1];
        pool->uz[k] = ev->neutron_dir[2];
        pool->e[k] = ev->neutron_e_mev;
        /* ev->local_deposit_mev: heavy recoil energy — deposit when scoring wired.
         * secondaries[0]: H-1 recoil proton — ion-feedback path not yet wired. */
        break;

    case OSH_NEUTRON_REACTION_COMPOUND:
        /* Push neutron secondaries (from FBU / (n,2n)) back for next pass.
         * Ion secondaries deposit locally until ion-feedback flag is wired. */
        for (i = 0u; i < ev->n_secondaries; ++i) {
            sec = &ev->secondaries[i];
            if (sec->species != NULL && sec->species->pdg == OSH_PART_PDG_NEUTRON) {
                push_neutron_secondary(pool, k, sec);
            }
        }
        pool->e[k] = 0.0;
        break;

    case OSH_NEUTRON_REACTION_CAPTURE:
    case OSH_NEUTRON_REACTION_CHARGE_EXCHANGE:
    case OSH_NEUTRON_REACTION_LOCAL_DEPOSIT:
    case OSH_NEUTRON_REACTION_NONE:
    default:
        /* Deposit ev->local_deposit_mev when scoring is wired.
         * CHARGE_EXCHANGE: secondary ion goes to ion pool with ion-feedback. */
        pool->e[k] = 0.0;
        break;
    }
}

/* --------------------------------------------------------------------------
 * Transport loop
 * -------------------------------------------------------------------------- */

enum osh_status osh_transport_neutron_run(struct osh_transport_context *transport_ctx,
                                          struct osh_beam_runtime *beam_rt,
                                          struct osh_gemca_runtime const *geom_rt,
                                          struct osh_material_runtime const *material_rt,
                                          struct osh_scoring_runtime *score_rt) {
    struct osh_neutron_pool *pool;             /* borrowed from transport_ctx */
    struct osh_nuclear_handler const *handler; /* element composition + FBU */
    struct osh_neutron_xsec xsec;              /* cross-section model (stack-alloc) */
    struct osh_zone_ref *zone_refs;            /* per-wavefront geometry scratch */
    double *dist_batch;                        /* per-wavefront boundary distance */
    size_t n_wavefront;                        /* snapshot of pool->n at pass start */
    size_t k;                                  /* slot index within wavefront */
    size_t mat_idx;                            /* material runtime index for slot k */
    double rho;                                /* material density [g/cm³] */
    double e_mev;                              /* neutron kinetic energy [MeV] */
    double e_cutoff_mev;                       /* effective neutron cutoff [MeV] */
    double dir[3];                             /* neutron direction unit vector */
    double sigma_tot;                          /* macroscopic total XS [cm⁻¹] */
    double l;                                  /* sampled free path length [cm] */
    struct osh_neutron_reaction_event ev;      /* reaction outcome */
    enum osh_status rc;
    size_t n_cutoff;
    size_t n_escape;
    size_t n_steps;
    size_t n_interactions;

    (void) beam_rt; /* no primary beam refill: neutrons come from the pool */
    if (!transport_ctx || !geom_rt || !material_rt) {
        return OSH_EINVAL;
    }

    pool = transport_ctx->neutron_pool;
    handler = transport_ctx->nuclear_handler;

    if (!pool || pool->n == 0u) {
        return OSH_OK; /* neutron pool absent or empty: nothing to do */
    }
    if (!handler) {
        OSH_DIAG_ERRORF(transport_ctx->diag, "%s", "neutron transport: nuclear_handler is required but not set");
        return OSH_EINVAL;
    }

    rc = osh_neutron_xsec_compile(transport_ctx->diag, &xsec);
    if (rc != OSH_OK) {
        return rc;
    }
    e_cutoff_mev = (transport_ctx->params.ncut > 0.0f) ? (double) transport_ctx->params.ncut
                                                       : (double) OSH_TRANSPORT_NEUTRON_CUTOFF_DEFAULT_MEV;
    n_cutoff = 0u;
    n_escape = 0u;
    n_steps = 0u;
    n_interactions = 0u;

    /* Borrow pre-allocated geometry scratch from the transport context.
     * Capacity was fixed at simulation create time; both pools share this scratch
     * since ion and neutron passes run sequentially. */
    zone_refs = transport_ctx->zone_refs;
    dist_batch = transport_ctx->dist_batch;
    if (!zone_refs || !dist_batch) {
        osh_neutron_xsec_free(&xsec);
        return OSH_EINVAL;
    }

    OSH_DIAG_INFOF(transport_ctx->diag, "Neutron transport: starting with %zu neutrons", pool->n);

    /* ---- Wavefront loop -------------------------------------------------- */
    while (pool->n > 0u) {
        /* Snapshot live count: secondaries pushed this pass go to slots
         * [n_wavefront..pool->n) and are processed in the next iteration. */
        n_wavefront = pool->n;

        osh_gemca_runtime_get_zone_ref_batch(
            geom_rt, pool->x, pool->y, pool->z, pool->ux, pool->uy, pool->uz, n_wavefront, zone_refs);
        osh_gemca_runtime_get_distance_batch(
            geom_rt, pool->x, pool->y, pool->z, pool->ux, pool->uy, pool->uz, zone_refs, n_wavefront, dist_batch);

        for (k = 0u; k < n_wavefront; ++k) {
            /* -- energy cutoff --------------------------------------------- */
            if (pool->e[k] <= e_cutoff_mev) {
                pool->e[k] = 0.0;
                n_cutoff++;
                continue;
            }

            /* -- escaped geometry ------------------------------------------ */
            if (zone_refs[k].zone_idx == OSH_ZONE_INDEX_INVALID) {
                pool->e[k] = 0.0;
                n_escape++;
                continue;
            }

            mat_idx = zone_refs[k].material_idx;

            /* -- blackhole -------------------------------------------------- */
            if (mat_idx == (size_t) OSH_MATERIAL_INDEX_BLACKHOLE) {
                pool->e[k] = 0.0;
                continue;
            }

            rho = (double) material_rt->rho[mat_idx];

            /* -- vacuum or zero-density material --------------------------- */
            if (rho <= 0.0) {
                rc = score_neutron_step(score_rt, pool, k, &zone_refs[k], rho, dist_batch[k]);
                if (rc != OSH_OK) {
                    osh_neutron_xsec_free(&xsec);
                    return rc;
                }
                n_steps++;
                advance_to_boundary(pool, k, dist_batch[k]);
                continue;
            }

            e_mev = pool->e[k];
            sigma_tot = neutron_sigma_tot_cm(&xsec, handler, mat_idx, rho, e_mev);

            /* -- transparent material (no xsec for any element) ------------ */
            if (sigma_tot <= 0.0) {
                rc = score_neutron_step(score_rt, pool, k, &zone_refs[k], rho, dist_batch[k]);
                if (rc != OSH_OK) {
                    osh_neutron_xsec_free(&xsec);
                    return rc;
                }
                n_steps++;
                advance_to_boundary(pool, k, dist_batch[k]);
                continue;
            }

            /* -- sample free path length ------------------------------------ */
            /* u=0 is possible (1-in-2^32 for PCG32) and log(0)=-inf raises
             * FE_DIVBYZERO on trapping platforms; HUGE_VAL gives inf/sigma_tot
             * = inf which correctly falls through to the boundary-crossing path. */
            {
                double u = osh_rng_double(&pool->rng[k]);
                l = (u > 0.0) ? (-log(u) / sigma_tot) : HUGE_VAL;
            }

            if (l >= dist_batch[k]) {
                /* no interaction this step: neutron crosses zone boundary */
                rc = score_neutron_step(score_rt, pool, k, &zone_refs[k], rho, dist_batch[k]);
                if (rc != OSH_OK) {
                    osh_neutron_xsec_free(&xsec);
                    return rc;
                }
                n_steps++;
                advance_to_boundary(pool, k, dist_batch[k]);
                continue;
            }

            /* -- advance to interaction point ------------------------------- */
            rc = score_neutron_step(score_rt, pool, k, &zone_refs[k], rho, l);
            if (rc != OSH_OK) {
                osh_neutron_xsec_free(&xsec);
                return rc;
            }
            n_steps++;
            n_interactions++;
            pool->x[k] += pool->ux[k] * l;
            pool->y[k] += pool->uy[k] * l;
            pool->z[k] += pool->uz[k] * l;

            dir[0] = pool->ux[k];
            dir[1] = pool->uy[k];
            dir[2] = pool->uz[k];

            /* -- sample and apply reaction ---------------------------------- */
            osh_neutron_reaction_sample(&xsec, handler, mat_idx, rho, e_mev, dir, &pool->rng[k], &ev);
            apply_event(pool, k, &ev);
        }

        osh_neutron_pool_compact(pool);
    }

    OSH_DIAG_INFOF(transport_ctx->diag,
                   "Neutron transport: %zu created, %zu dropped, %zu steps, %zu interactions, %zu cutoff, %zu escaped",
                   pool->n_created,
                   pool->n_dropped,
                   n_steps,
                   n_interactions,
                   n_cutoff,
                   n_escape);

    osh_neutron_xsec_free(&xsec);
    return OSH_OK;
}
