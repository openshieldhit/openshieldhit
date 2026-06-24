#include "physics/neutron/osh_neutron_reaction.h"

#include <math.h>
#include <string.h>

#include "openshieldhit/const.h"
#include "particle/osh_particle.h"
#include "particle/osh_particle_const.h"
#include "particle/osh_particle_pdg.h"
#include "physics/neutron/osh_neutron_xsec.h"
#include "physics/nuclear/osh_nuclear_compound.h"
#include "physics/osh_kinematics.h"
#include "random/osh_rng.h"

/* Maximum number of nuclides per material supported in this layer. */
#define MAX_ELEMS 64

/* --------------------------------------------------------------------------
 * Species cache (proton, alpha) for CHARGE_EXCHANGE secondaries.
 * Filled once on first call; all transport is single-threaded.
 * -------------------------------------------------------------------------- */

static struct particle s_proton;
static struct particle s_alpha;
static int s_species_ready = 0;

static void ensure_species(void) {
    if (s_species_ready)
        return;
    osh_particle_from_pdg(&s_proton, OSH_PART_PDG_PROTON);
    osh_particle_from_pdg(&s_alpha, OSH_PART_PDG_HE4);
    s_species_ready = 1;
}

/* Fill event with a local energy deposit and no secondaries. */
static void set_local_deposit(struct osh_neutron_reaction_event *ev, double e_mev) {
    ev->kind = OSH_NEUTRON_REACTION_LOCAL_DEPOSIT;
    ev->local_deposit_mev = e_mev;
    ev->n_secondaries = 0u;
}

/* --------------------------------------------------------------------------
 * Elastic scatter kinematics (non-relativistic, neutron mass ≈ target mass/A)
 * -------------------------------------------------------------------------- */

static void do_elastic(
    unsigned int a, double e_mev, double const dir[3], struct osh_rng *rng, struct osh_neutron_reaction_event *ev) {
    double cos_theta_cm;  /* CM scatter angle cosine (sampled isotropically — P0 approx) */
    double cos_theta_lab; /* lab neutron scatter angle cosine, from NR frame transform */
    double sin_theta_lab; /* derived; passed to rotate_dir_cos */
    double cos_phi;       /* azimuthal direction (no trig; Knuth disk-rejection) */
    double sin_phi;
    double af;    /* A as double — appears repeatedly in NR energy formula */
    double denom; /* (1 + A)², common denominator in NR kinematics */

    af = (double) a;

    /* isotropic CM (P0 — no JEFF angular distribution for minimal model) */
    cos_theta_cm = 2.0 * osh_rng_double(rng) - 1.0;
    osh_kinematics_azimuth(rng, &cos_phi, &sin_phi);

    /* NR scattered neutron energy: E_n' = T_n (A² + 1 + 2A cosθ_CM) / (A+1)² */
    denom = (1.0 + af) * (1.0 + af);
    ev->neutron_e_mev = e_mev * (af * af + 1.0 + 2.0 * af * cos_theta_cm) / denom;
    if (ev->neutron_e_mev < 0.0)
        ev->neutron_e_mev = 0.0;

    /* NR lab scatter angle: cos θ_lab = (1 + A cosθ_CM) / sqrt(A² + 2A cosθ_CM + 1)
     * Guard: denominator is zero when A=1 and cosθ_CM = -1 (exact backward CM scatter),
     * which also gives E_n' = 0.  Direction is meaningless at zero energy; copy the
     * incident direction and let the energy cutoff remove the neutron. */
    if (ev->neutron_e_mev > 0.0) {
        cos_theta_lab = (1.0 + af * cos_theta_cm) / sqrt(af * af + 2.0 * af * cos_theta_cm + 1.0);
        sin_theta_lab = sqrt(fmax(0.0, 1.0 - cos_theta_lab * cos_theta_lab));
        osh_kinematics_rotate_dir_cos(dir, ev->neutron_dir, cos_theta_lab, sin_theta_lab, cos_phi, sin_phi);
    } else {
        ev->neutron_dir[0] = dir[0];
        ev->neutron_dir[1] = dir[1];
        ev->neutron_dir[2] = dir[2];
    }

    ev->kind = OSH_NEUTRON_REACTION_ELASTIC;
    ev->n_secondaries = 0u;
    /* heavy recoil below transport threshold: deposit locally.
     * Future: push to fragment pool and score with score_point(). */
    ev->local_deposit_mev = fmax(0.0, e_mev - ev->neutron_e_mev);

    /* H-1 recoil proton: up to full energy transfer, important for dose in tissue */
    if (a == 1u) {
        double t_recoil;   /* proton kinetic energy = T_n - T_n' (NR energy conservation) */
        double cos_recoil; /* NR relation: cos θ_recoil = sqrt((1 - cosθ_CM) / 2) */
        double sin_recoil;
        struct osh_nuclear_secondary *sec;

        t_recoil = e_mev - ev->neutron_e_mev;
        if (t_recoil > 0.0) {
            cos_recoil = sqrt(fmax(0.0, (1.0 - cos_theta_cm) / 2.0));
            sin_recoil = sqrt(fmax(0.0, 1.0 - cos_recoil * cos_recoil));
            sec = &ev->secondaries[0];
            /* recoil azimuth is exactly opposite to scattered neutron */
            osh_kinematics_rotate_dir_cos(dir, sec->dir, cos_recoil, sin_recoil, -cos_phi, -sin_phi);
            sec->energy = t_recoil;
            sec->species = &s_proton;
            ev->n_secondaries = 1u;
            ev->local_deposit_mev = 0.0; /* proton secondary carries the recoil */
        }
    }
}

/* --------------------------------------------------------------------------
 * 2-body charge-exchange kinematics: n + (Z,A) → product + residual
 * Full relativistic via osh_kinematics_boost_to_lab.
 * -------------------------------------------------------------------------- */

static void do_charge_exchange(unsigned int z,
                               unsigned int a,
                               int product_pdg,
                               unsigned int z_res,
                               unsigned int a_res,
                               double e_mev,
                               double const dir[3],
                               struct osh_rng *rng,
                               struct osh_neutron_reaction_event *ev) {
    double m_target;     /* nuclear rest mass M(Z,A) [MeV/c²] */
    double m_product;    /* rest mass of outgoing light ion (p or α) [MeV/c²] */
    double m_residual;   /* rest mass of recoil nucleus M(Z_res, A_res) [MeV/c²] */
    double p_n_mag;      /* |p_n| = sqrt(T² + 2Tm_n) [MeV/c], relativistic */
    double p_parent[3];  /* lab momentum of (n + target) system = p_n (target at rest) */
    double W;            /* invariant mass sqrt(s) of n+target system [MeV/c²] */
    double p_cm_mag;     /* CM momentum of each product from Källén triangle */
    double cos_theta_cm; /* isotropic CM direction for the light product */
    double sin_theta_cm;
    double cos_phi;
    double sin_phi;
    double p_cm[3];   /* light product momentum in CM frame */
    double E1_cm;     /* total energy of light product in CM frame [MeV] */
    double e_lab;     /* total energy of light product in lab [MeV] */
    double p_lab[3];  /* lab momentum of light product [MeV/c] */
    double p_lab_mag; /* magnitude of p_lab, for direction normalisation */
    double T_product; /* kinetic energy of outgoing ion in lab [MeV] */
    struct osh_nuclear_secondary *sec;

    if (!osh_particle_nuclear_mass_mev_from_za(z, a, &m_target)) {
        set_local_deposit(ev, e_mev);
        return;
    }
    if (!osh_particle_nuclear_mass_mev_from_za(z_res, a_res, &m_residual)) {
        set_local_deposit(ev, e_mev);
        return;
    }
    if (!osh_particle_nuclear_mass_from_pdg(product_pdg, &m_product)) {
        set_local_deposit(ev, e_mev);
        return;
    }

    /* lab momentum of incident neutron (along dir) */
    p_n_mag = sqrt(e_mev * e_mev + 2.0 * e_mev * OSH_PART_MASS_NEUTRON);
    p_parent[0] = p_n_mag * dir[0];
    p_parent[1] = p_n_mag * dir[1];
    p_parent[2] = p_n_mag * dir[2];

    /* invariant mass of (n + target) system */
    W = sqrt((OSH_PART_MASS_NEUTRON + m_target) * (OSH_PART_MASS_NEUTRON + m_target) + 2.0 * e_mev * m_target);

    /* CM momentum of products */
    p_cm_mag = osh_kinematics_two_body_decay_p(W, m_product, m_residual);
    if (p_cm_mag <= 0.0) {
        /* channel closed (Q < 0 and T too low) */
        set_local_deposit(ev, e_mev);
        return;
    }

    /* isotropic CM direction */
    cos_theta_cm = 2.0 * osh_rng_double(rng) - 1.0;
    sin_theta_cm = sqrt(fmax(0.0, 1.0 - cos_theta_cm * cos_theta_cm));
    osh_kinematics_azimuth(rng, &cos_phi, &sin_phi);
    p_cm[0] = p_cm_mag * sin_theta_cm * cos_phi;
    p_cm[1] = p_cm_mag * sin_theta_cm * sin_phi;
    p_cm[2] = p_cm_mag * cos_theta_cm;

    /* product total energy in CM */
    E1_cm = sqrt(m_product * m_product + p_cm_mag * p_cm_mag);

    /* boost to lab */
    osh_kinematics_boost_to_lab(W, p_parent, E1_cm, p_cm, &e_lab, p_lab);
    T_product = e_lab - m_product;
    if (T_product < 0.0)
        T_product = 0.0;

    p_lab_mag = sqrt(p_lab[0] * p_lab[0] + p_lab[1] * p_lab[1] + p_lab[2] * p_lab[2]);

    sec = &ev->secondaries[0];
    if (p_lab_mag > 0.0) {
        sec->dir[0] = p_lab[0] / p_lab_mag;
        sec->dir[1] = p_lab[1] / p_lab_mag;
        sec->dir[2] = p_lab[2] / p_lab_mag;
    } else {
        sec->dir[0] = dir[0];
        sec->dir[1] = dir[1];
        sec->dir[2] = dir[2];
    }
    sec->energy = T_product;
    sec->species = (product_pdg == OSH_PART_PDG_PROTON) ? &s_proton : &s_alpha;

    ev->kind = OSH_NEUTRON_REACTION_CHARGE_EXCHANGE;
    ev->n_secondaries = 1u;
    ev->local_deposit_mev = 0.0;
}

/* --------------------------------------------------------------------------
 * Compound nucleus route: build (Z, A+1, E*) and dispatch
 * -------------------------------------------------------------------------- */

static void do_compound(unsigned int z,
                        unsigned int a,
                        double e_mev,
                        double const dir[3],
                        struct osh_nuclear_handler const *handler,
                        struct osh_diag_sink const *diag,
                        struct osh_rng *rng,
                        struct osh_neutron_reaction_event *ev) {
    double m_target;                 /* nuclear rest mass M(Z,A) [MeV/c²] */
    double m_compound;               /* nuclear rest mass M(Z,A+1) [MeV/c²] */
    double s_n;                      /* neutron separation energy: S_n = m_n + M(Z,A) - M(Z,A+1) */
    double e_star;                   /* excitation energy: E* = S_n + T_n · A/(A+1) */
    double p_n_mag;                  /* incident neutron lab momentum magnitude [MeV/c] */
    double p_lab[3];                 /* compound nucleus lab momentum = p_n (target at rest) */
    unsigned int z_c;                /* compound nucleus atomic number (= target Z) */
    unsigned int a_c;                /* compound nucleus mass number (= target A + 1) */
    struct osh_nuclear_event fbu_ev; /* output from FBU or heavy-A sink */
    size_t i;

    z_c = z;
    a_c = a + 1u;

    /* neutron separation energy S_n = m_n + M(Z,A) - M(Z,A+1) */
    if (!osh_particle_nuclear_mass_mev_from_za(z, a, &m_target)
        || !osh_particle_nuclear_mass_mev_from_za(z_c, a_c, &m_compound)) {
        set_local_deposit(ev, e_mev);
        return;
    }
    s_n = OSH_PART_MASS_NEUTRON + m_target - m_compound;

    /* excitation energy: E* = S_n + T_n * A/(A+1) */
    e_star = s_n + e_mev * (double) a / (double) a_c;
    if (e_star < 0.0)
        e_star = 0.0;

    /* compound nucleus lab momentum = incident neutron momentum */
    p_n_mag = sqrt(e_mev * e_mev + 2.0 * e_mev * OSH_PART_MASS_NEUTRON);
    p_lab[0] = p_n_mag * dir[0];
    p_lab[1] = p_n_mag * dir[1];
    p_lab[2] = p_n_mag * dir[2];

    osh_nuclear_compound_step(z_c, a_c, e_star, p_lab, &handler->fbu, diag, rng, &fbu_ev);

    /* translate FBU/sink event into neutron reaction event */
    if (fbu_ev.kind == OSH_NUCLEAR_EVENT_ABSORB || fbu_ev.n_secondaries == 0u) {
        /* heavy-A sink or FBU produced nothing transportable */
        set_local_deposit(ev, e_star);
        return;
    }

    ev->kind = OSH_NEUTRON_REACTION_COMPOUND;
    ev->local_deposit_mev = 0.0;
    ev->n_secondaries =
        fbu_ev.n_secondaries > OSH_NUCLEAR_MAX_SECONDARIES ? OSH_NUCLEAR_MAX_SECONDARIES : fbu_ev.n_secondaries;

    for (i = 0u; i < ev->n_secondaries; ++i)
        ev->secondaries[i] = fbu_ev.secondaries[i];
}

/* --------------------------------------------------------------------------
 * Public API
 * -------------------------------------------------------------------------- */

void osh_neutron_reaction_sample(struct osh_neutron_xsec *xsec,
                                 struct osh_nuclear_handler const *handler,
                                 size_t material_idx,
                                 double rho_g_cm3,
                                 double e_mev,
                                 double const dir[3],
                                 struct osh_rng *rng,
                                 struct osh_neutron_reaction_event *event_out) {
    struct osh_nuclear_elem const *elems; /* element list for this material */
    size_t n_elems;
    struct osh_neutron_xsec_result sig[MAX_ELEMS]; /* per-element cross-section cache */
    double nd_times_tot[MAX_ELEMS];                /* nᵢ · σ_tot,i [mb·cm⁻³]: weight for target sampling */
    double cum_tot;                                /* sum of nd_times_tot; normalises the target roulette */
    double xi;                                     /* uniform deviate in [0, cum_tot) for target element selection */
    double running;                                /* partial sum during target element scan */
    size_t i;
    size_t chosen; /* index of sampled target element */
    unsigned int z;
    unsigned int a;
    double xi2; /* uniform deviate in [0, σ_tot) for reaction channel selection */
    double cum; /* cumulative channel cross section during channel roulette */

    memset(event_out, 0, sizeof(*event_out));
    ensure_species();

    /* -- 1. element list for this material --------------------------------- */
    elems = handler->elem_pool + handler->elem_offset[material_idx];
    n_elems = handler->elem_count[material_idx];
    if (n_elems == 0u || n_elems > MAX_ELEMS) {
        event_out->kind = OSH_NEUTRON_REACTION_NONE;
        return;
    }

    /* -- 2. per-element number density × σ_tot (mb·cm⁻³) ------------------ */
    cum_tot = 0.0;
    for (i = 0u; i < n_elems; ++i) {
        double nd_i;
        /* n_i [cm⁻³] = w_i * ρ [g/cm³] * N_A [mol⁻¹] / A_i [g/mol] */
        nd_i = (double) elems[i].mass_fraction * rho_g_cm3 * OSH_NAVOGADRO / (double) elems[i].a;
        osh_neutron_xsec_lookup(xsec, (int) elems[i].z, (int) elems[i].a, e_mev, &sig[i]);
        nd_times_tot[i] = nd_i * sig[i].tot; /* mb·cm⁻³ */
        cum_tot += nd_times_tot[i];
    }

    if (cum_tot <= 0.0) {
        event_out->kind = OSH_NEUTRON_REACTION_NONE;
        return;
    }

    /* -- 3. sample target proportional to nᵢ σ_tot,i ---------------------- */
    xi = osh_rng_double(rng) * cum_tot;
    chosen = n_elems - 1u; /* fallback: last element */
    running = 0.0;
    for (i = 0u; i < n_elems; ++i) {
        running += nd_times_tot[i];
        if (running >= xi) {
            chosen = i;
            break;
        }
    }

    z = elems[chosen].z;
    a = elems[chosen].a;

    /* -- 4. sample reaction channel ---------------------------------------- */
    xi2 = osh_rng_double(rng) * sig[chosen].tot;
    cum = 0.0;

    /* elastic */
    cum += sig[chosen].el;
    if (xi2 < cum) {
        do_elastic(a, e_mev, dir, rng, event_out);
        return;
    }

    /* (n,γ) capture */
    cum += sig[chosen].ng;
    if (xi2 < cum) {
        event_out->kind = OSH_NEUTRON_REACTION_CAPTURE;
        event_out->local_deposit_mev = e_mev;
        event_out->n_secondaries = 0u;
        return;
    }

    /* (n,p): residual is (Z-1, A) */
    cum += sig[chosen].np;
    if (xi2 < cum && z >= 1u) {
        do_charge_exchange(z, a, OSH_PART_PDG_PROTON, z - 1u, a, e_mev, dir, rng, event_out);
        return;
    }

    /* (n,α): residual is (Z-2, A-3) */
    cum += sig[chosen].na;
    if (xi2 < cum && z >= 2u && a >= 4u) {
        do_charge_exchange(z, a, OSH_PART_PDG_HE4, z - 2u, a - 3u, e_mev, dir, rng, event_out);
        return;
    }

    /* (n,n') + (n,2n) + remainder → compound nucleus */
    do_compound(z, a, e_mev, dir, handler, xsec->diag, rng, event_out);
}
