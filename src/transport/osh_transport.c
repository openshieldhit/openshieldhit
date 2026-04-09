#include "transport/osh_transport.h"

#include "beam/osh_beam.h"
#include "beam/osh_beam_model.h"
#include "beam/osh_beamdef.h"
#include "common/osh_coord.h"
#include "gemca/osh_gemca2.h"
#include "material/osh_material.h"
#include "material/runtime/osh_material_runtime.h"
#include "random/osh_rng.h"
#include "scoring/runtime/osh_scoring_step.h"

#define OSH_TRANSPORT_BOUNDARY_EPS 1e-8
#define OSH_TRANSPORT_MAX_STEPS_PER_PRIMARY 1000000u

static enum osh_status transport_one_primary(struct gemca_workspace *geom,
                                             struct beam_workspace const *beam,
                                             struct material_workspace const *materials,
                                             struct osh_material_runtime const *tables,
                                             struct osh_scoring_runtime *scoring,
                                             struct particle const *part,
                                             struct ray_v *state,
                                             double deltae);
static double cutoff_total_energy(struct beam_workspace const *beam,
                                  struct osh_material_runtime const *tables,
                                  struct particle const *part);
static double energy_from_residual_range(struct osh_material_runtime const *tables,
                                         size_t material_idx,
                                         size_t projectile_idx,
                                         double residual_range);
static double energy_grid_value(struct osh_material_runtime const *tables, size_t energy_idx);
static enum osh_status find_projectile_index(struct osh_material_runtime const *tables,
                                             struct particle const *part,
                                             size_t *projectile_idx_out);
static int is_blackhole_material(size_t material_idx);
static int is_vacuum_material(size_t material_idx);
static void ray_from_state(struct ray *ray, struct ray_v const *state);
static void step_from_state(
    struct step *st, struct ray_v const *state, double step_len, double exit_energy, double rho, int medium, int zone);
static void advance_state(struct ray_v *state, double step_len, double exit_energy);
static void nudge_state(struct ray_v *state, double eps);

enum osh_status osh_transport_run_minimal(struct beam_workspace const *beam,
                                          struct gemca_workspace *geom,
                                          struct material_workspace const *materials,
                                          struct osh_material_runtime const *tables,
                                          struct osh_scoring_runtime *scoring) {
    struct osh_rng rng;
    struct particle *part;
    struct ray_v state;
    size_t iprim;
    enum osh_status rc;

    if (!beam || !geom || !materials || !tables || !scoring) {
        return OSH_EINVAL;
    }
    if (beam->beam_mode != OSH_BEAM_MODE_SPOTS) {
        return OSH_ENOTSUP;
    }
    if (beam->nstat == 0u) {
        return OSH_EINVAL;
    }
    if (beam->deltae <= 0.0f || beam->deltae >= 1.0f) {
        return OSH_EINVAL;
    }

    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, (uint64_t) beam->rndseed, (uint64_t) beam->rndoffset);

    for (iprim = 0; iprim < beam->nstat; ++iprim) {
        part = NULL;
        rc = osh_beam_new_primary(beam, &rng, &part, &state);
        if (rc != OSH_OK) {
            return rc;
        }
        rc = transport_one_primary(geom, beam, materials, tables, scoring, part, &state, (double) beam->deltae);
        if (rc != OSH_OK) {
            return rc;
        }
    }

    return OSH_OK;
}

static enum osh_status transport_one_primary(struct gemca_workspace *geom,
                                             struct beam_workspace const *beam,
                                             struct material_workspace const *materials,
                                             struct osh_material_runtime const *tables,
                                             struct osh_scoring_runtime *scoring,
                                             struct particle const *part,
                                             struct ray_v *state,
                                             double deltae) {
    struct ray ray;
    size_t projectile_idx;
    double cutoff_total;
    unsigned int istep;
    enum osh_status rc;

    if (!geom || !materials || !tables || !scoring || !part || !state) {
        return OSH_EINVAL;
    }

    rc = find_projectile_index(tables, part, &projectile_idx);
    if (rc != OSH_OK) {
        return rc;
    }
    cutoff_total = cutoff_total_energy(beam, tables, part);

    for (istep = 0; istep < OSH_TRANSPORT_MAX_STEPS_PER_PRIMARY; ++istep) {
        size_t zone_idx;
        struct zone *zone;
        struct material const *material;
        double boundary_ds;
        double rho;
        double step_len;
        double exit_energy;
        struct step st;
        int hit_boundary;

        if (state->p[3] <= cutoff_total) {
            return OSH_OK;
        }

        ray_from_state(&ray, state);
        zone_idx = osh_gemca_get_zone_index(geom, &ray);
        if (zone_idx == OSH_GEMCA_ZONE_INDEX_INVALID) {
            return OSH_OK;
        }

        zone = geom->zones[zone_idx];
        if (!zone) {
            return OSH_ESTATE;
        }
        if (is_blackhole_material(zone->material_idx)) {
            return OSH_OK;
        }

        material = osh_material_by_index(materials, zone->material_idx);
        if (!material) {
            return OSH_ESTATE;
        }

        boundary_ds = osh_gemca_get_distance(zone, &ray);
        if (boundary_ds < 0.0) {
            return OSH_ESTATE;
        }
        if (boundary_ds <= OSH_TRANSPORT_BOUNDARY_EPS) {
            nudge_state(state, OSH_TRANSPORT_BOUNDARY_EPS);
            continue;
        }

        rho = (material->rho > 0.0) ? material->rho : 0.0;
        if (is_vacuum_material(zone->material_idx) || rho <= 0.0) {
            step_len = boundary_ds;
            exit_energy = state->p[3];
        } else {
            double a_proj;
            double e0_total;
            double e0_per_nuc;
            double e1_target;
            double r0;
            double r1_target;
            double ds_csda;
            double residual_range;

            a_proj = (part->a > 0u) ? (double) part->a : 1.0;
            e0_total = state->p[3];
            e0_per_nuc = e0_total / a_proj;
            e1_target = e0_total * (1.0 - deltae);
            if (e1_target < cutoff_total) {
                e1_target = cutoff_total;
            }
            if (e1_target >= e0_total) {
                return OSH_OK;
            }

            r0 = osh_material_runtime_range_lookup(tables, zone->material_idx, projectile_idx, e0_per_nuc);
            r1_target =
                osh_material_runtime_range_lookup(tables, zone->material_idx, projectile_idx, e1_target / a_proj);
            ds_csda = (r0 - r1_target) / rho;
            if (ds_csda <= 0.0) {
                return OSH_OK;
            }

            if (ds_csda <= boundary_ds) {
                step_len = ds_csda;
                exit_energy = e1_target;
            } else {
                step_len = boundary_ds;
                residual_range = r0 - rho * step_len;
                if (residual_range < 0.0) {
                    residual_range = 0.0;
                }
                exit_energy =
                    energy_from_residual_range(tables, zone->material_idx, projectile_idx, residual_range) * a_proj;
                if (exit_energy < cutoff_total) {
                    exit_energy = cutoff_total;
                }
            }
        }

        if (step_len <= 0.0) {
            return OSH_OK;
        }

        if (exit_energy > state->p[3]) {
            exit_energy = state->p[3];
        }

        step_from_state(&st, state, step_len, exit_energy, rho, (int) zone->material_idx, (int) zone_idx);
        rc = osh_scoring_score_step(scoring, part, &st);
        if (rc != OSH_OK) {
            return rc;
        }

        hit_boundary = (step_len >= boundary_ds - OSH_TRANSPORT_BOUNDARY_EPS);
        advance_state(state, step_len, exit_energy);
        if (hit_boundary) {
            nudge_state(state, OSH_TRANSPORT_BOUNDARY_EPS);
        }
    }

    return OSH_ESTATE;
}

static double cutoff_total_energy(struct beam_workspace const *beam,
                                  struct osh_material_runtime const *tables,
                                  struct particle const *part) {
    double a_proj;
    double cutoff_total;
    double cutoff_from_transport;
    double cutoff_from_beam;

    a_proj = (part->a > 0u) ? (double) part->a : 1.0;
    cutoff_total = OSH_BEAM_TMIN;
    cutoff_from_beam = 0.0;
    if (beam) {
        cutoff_from_beam = (double) beam->tcut * a_proj;
        if (cutoff_from_beam > cutoff_total) {
            cutoff_total = cutoff_from_beam;
        }
    }
    cutoff_from_transport = tables->emin * a_proj;
    if (cutoff_from_transport > cutoff_total) {
        cutoff_total = cutoff_from_transport;
    }
    return cutoff_total;
}

static double energy_from_residual_range(struct osh_material_runtime const *tables,
                                         size_t material_idx,
                                         size_t projectile_idx,
                                         double residual_range) {
    float const *range_col;
    size_t lo;
    size_t hi;
    size_t mid;
    double r_lo;
    double r_hi;
    double e_lo;
    double e_hi;
    double frac;

    range_col = tables->range_csda + (material_idx * tables->nprojectiles + projectile_idx) * tables->nenergy;
    if (residual_range <= (double) range_col[0]) {
        return tables->emin;
    }
    if (residual_range >= (double) range_col[tables->nenergy - 1u]) {
        return tables->emax;
    }

    lo = 0u;
    hi = tables->nenergy - 1u;
    while (hi - lo > 1u) {
        mid = lo + (hi - lo) / 2u;
        if ((double) range_col[mid] <= residual_range) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    r_lo = (double) range_col[lo];
    r_hi = (double) range_col[hi];
    e_lo = energy_grid_value(tables, lo);
    e_hi = energy_grid_value(tables, hi);
    if (r_hi <= r_lo) {
        return e_lo;
    }
    frac = (residual_range - r_lo) / (r_hi - r_lo);
    return e_lo * (1.0 - frac) + e_hi * frac;
}

static double energy_grid_value(struct osh_material_runtime const *tables, size_t energy_idx) {
    return exp(tables->log_emin + (double) energy_idx / tables->inv_dlog);
}

static enum osh_status find_projectile_index(struct osh_material_runtime const *tables,
                                             struct particle const *part,
                                             size_t *projectile_idx_out) {
    unsigned int z_match;
    unsigned int a_match;
    size_t i;

    if (!tables || !part || !projectile_idx_out) {
        return OSH_EINVAL;
    }
    z_match = part->z;
    a_match = part->a;
    if (part->pdg == 2212 && z_match == 0u && a_match == 0u) {
        z_match = 1u;
        a_match = 1u;
    }
    for (i = 0; i < tables->nprojectiles; ++i) {
        if (tables->projectile_z[i] == z_match && tables->projectile_a[i] == a_match) {
            *projectile_idx_out = i;
            return OSH_OK;
        }
    }
    for (i = 0; i < tables->nprojectiles; ++i) {
        if (tables->projectile_z[i] == z_match) {
            *projectile_idx_out = i;
            return OSH_OK;
        }
    }
    return OSH_ENOTSUP;
}

static int is_blackhole_material(size_t material_idx) {
    return material_idx == OSH_MATERIAL_INDEX_BLACKHOLE;
}

static int is_vacuum_material(size_t material_idx) {
    return material_idx == OSH_MATERIAL_INDEX_VACUUM;
}

static void ray_from_state(struct ray *ray, struct ray_v const *state) {
    ray->p[0] = state->p[0];
    ray->p[1] = state->p[1];
    ray->p[2] = state->p[2];
    ray->cp[0] = state->v[0];
    ray->cp[1] = state->v[1];
    ray->cp[2] = state->v[2];
    ray->system = state->system;
}

static void step_from_state(
    struct step *st, struct ray_v const *state, double step_len, double exit_energy, double rho, int medium, int zone) {
    int i;

    for (i = 0; i < 3; ++i) {
        st->p[i] = state->p[i];
        st->q[i] = state->p[i] + state->v[i] * step_len;
        st->v[i] = state->v[i];
        st->w[i] = state->v[i];
    }
    st->p[3] = state->p[3];
    st->q[3] = exit_energy;
    st->ds = step_len;
    st->de = state->p[3] - exit_energy;
    if (st->de < 0.0) {
        st->de = 0.0;
    }
    st->rho = rho;
    st->medium = medium;
    st->zone = zone;
    st->system = state->system;
}

static void advance_state(struct ray_v *state, double step_len, double exit_energy) {
    int i;

    for (i = 0; i < 3; ++i) {
        state->p[i] += state->v[i] * step_len;
    }
    state->p[3] = exit_energy;
}

static void nudge_state(struct ray_v *state, double eps) {
    int i;

    for (i = 0; i < 3; ++i) {
        state->p[i] += state->v[i] * eps;
    }
}
