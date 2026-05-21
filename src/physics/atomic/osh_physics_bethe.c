#include "physics/atomic/osh_physics_bethe.h"

#include <math.h>

#include "common/osh_const.h"

/* ---- Constants ----------------------------------------------------------- */

/*
 * K = 4 pi N_A r_e^2 m_e c^2 = 0.307075 MeV cm^2 / mol.
 *
 * This is the prefactor of the Bethe-Bloch formula.  It is exact to the
 * precision of the fundamental constants used by both libdedx (and SHIELD-HIT).
 * Reference: PDG Review of Particle Physics, "Passage of Particles through
 * Matter", Table 34.1.
 */
#define BETHE_K 0.307075

/*
 * 2 m_e c^2 = 1.022e6 eV.
 *
 * Appears in the argument of the Bethe logarithm as the numerator of
 * (2 m_e c^2 beta^2 gamma^2 / I), where I is in eV.
 * Value: 2 * 0.511 MeV = 1.022 MeV = 1.022e6 eV.
 */
#define BETHE_2MEC2_EV 1.022e6

/*
 * *** ENERGY UNIT AND MASS CONVENTION — READ CAREFULLY ***
 *
 * The public interface (osh_physics_bethe_eval, the sewing-point search, and
 * the LOADDEDX table axis) all use:
 *
 *   T/A  [MeV/nucleon]   where A is the INTEGER mass number of the ion.
 *
 * This is a counting convention: total kinetic energy T divided by the number
 * of nucleons A.  It is NOT the same as T/u (MeV per atomic mass unit),
 * although the two differ by less than 1% because nuclear binding energies
 * (~8 MeV/nucleon) are small compared to the nucleon rest mass (~931 MeV/c²).
 *
 * Internally, bethe_raw() receives total kinetic energy T = (T/A) * A [MeV]
 * and uses proj->mass_mev directly for the projectile rest mass.  Callers
 * must set proj->mass_mev to the NUCLEAR (fully-stripped) rest mass:
 *
 *   M_nuclear = amass [amu] * OSH_AMU  -  Z * m_electron
 *
 * This is provided by osh_particle_nuclear_mass_mev_from_za() in particle/.
 * Do NOT use A * 940 (free-nucleon approximation) or A * OSH_AMU (ignores
 * electron masses): both introduce ~1% error in the kinematic terms.
 *
 * For CSDA range: R [g/cm²] = ∫ dT / SP(T) = ∫ A · d(T/A) / SP(T/A).
 * The factor A must be included when integrating over the T/A axis; this
 * is handled by the a_proj argument to integrate_range() in prepare.c.
 */

/*
 * ln(10) = 4.60517... ~ 4.606.
 *
 * The Sternheimer density correction uses log10 for the momentum variable x
 * but natural log for the plasma frequency c.  The constant ln(10) converts
 * between the two.  The rounded value 4.606 is used throughout [SP].
 */
#define BETHE_LN10 4.606

/*
 * Golden-section search fractions.
 *
 * phi = (1 + sqrt(5)) / 2 ~ 1.6180339 (golden ratio).
 * The search bracket is split at the two interior points:
 *   lower fraction: 1 - 1/phi = (3 - sqrt(5)) / 2 ~ 0.3819661
 *   upper fraction:     1/phi = (sqrt(5) - 1) / 2 ~ 0.6180339
 * Using named constants avoids magic literals across the three search phases.
 */
#define GOLD_LOWER 0.3819661
#define GOLD_UPPER 0.6180339

/* ---- Forward declarations ------------------------------------------------ */

static double bethe_raw(double t_total_mev,
                        struct osh_physics_bethe_projectile const *proj,
                        struct osh_physics_bethe_target const *target);

/* ---- Public API ---------------------------------------------------------- */

void osh_physics_bethe_sewn_compute(struct osh_physics_bethe_projectile const *proj,
                                    struct osh_physics_bethe_target const *target,
                                    struct osh_physics_bethe_sewn *sewn) {
    /*
     * Golden-section search for the Lindhard-Scharff sewing point.
     *
     * The Bethe formula is only valid above a certain kinetic energy.  Below
     * that energy the stopping power is better described by the
     * Lindhard-Scharff (LSS) formula:
     *
     *   dE/dx_LSS(T) = C_lss * sqrt(T)
     *
     * where C_lss is chosen for continuity at the sewing point e_sewn.  The
     * sewing point is defined as the energy at which the tangent to the Bethe
     * curve passes through the origin when plotted as dE/dx vs sqrt(T), i.e.
     * the condition:
     *
     *   d/dT [dE/dx_Bethe(T)] = dE/dx_Bethe(T) / T
     *
     * equivalently:  d/dT [dE/dx / sqrt(T)] = 0.
     *
     * The search is identical to the three-phase golden-section procedure in
     * libdedx gold_section() [dedx_bethe.c].  Phase 1 finds the zero of the
     * Bethe formula (below which it becomes unphysical), phase 2 finds the
     * Bragg-peak maximum, and phase 3 finds the sewing point between those
     * two extremes.
     *
     * Search bounds are [1e-5, 10] MeV/nucleon.  Convergence tolerance 1e-5.
     */

    double e_min;
    double e_max;
    double e_zero;
    double e_extr;
    double rla;
    double rmu;
    double tla;
    double tmu;
    double epsilon;
    double h;
    double a;

    a = proj->a;
    epsilon = 1e-5;
    h = 1e-5;

    /* Phase 1: find e_zero, the energy per nucleon where Bethe crosses zero. */
    e_min = 1e-5;
    e_max = 1e1;
    rla = e_min + GOLD_LOWER * (e_max - e_min);
    rmu = e_min + GOLD_UPPER * (e_max - e_min);
    tla = bethe_raw(rla * a, proj, target);
    tmu = bethe_raw(rmu * a, proj, target);

    while ((e_max - e_min) >= epsilon) {
        if (tla > 0.0) {
            e_max = rmu;
            rmu = rla;
            tmu = tla;
            rla = e_min + GOLD_LOWER * (e_max - e_min);
            tla = bethe_raw(rla * a, proj, target);
        } else {
            e_min = rla;
            rla = rmu;
            tla = tmu;
            rmu = e_min + GOLD_UPPER * (e_max - e_min);
            tmu = bethe_raw(rmu * a, proj, target);
        }
    }
    e_zero = (e_max + e_min) * 0.5;

    /* Phase 2: find e_extr, the energy per nucleon of the Bragg-peak maximum,
     * searching [e_zero, 10] MeV/nucleon. */
    e_min = e_zero;
    e_max = 1e1;
    rla = e_min + GOLD_LOWER * (e_max - e_min);
    rmu = e_min + GOLD_UPPER * (e_max - e_min);
    tla = bethe_raw(rla * a, proj, target);
    tmu = bethe_raw(rmu * a, proj, target);

    while ((e_max - e_min) >= epsilon) {
        if (tla > tmu) {
            e_max = rmu;
            rmu = rla;
            tmu = tla;
            rla = e_min + GOLD_LOWER * (e_max - e_min);
            tla = bethe_raw(rla * a, proj, target);
        } else {
            e_min = rla;
            rla = rmu;
            tla = tmu;
            rmu = e_min + GOLD_UPPER * (e_max - e_min);
            tmu = bethe_raw(rmu * a, proj, target);
        }
    }
    e_extr = (e_min + e_max) * 0.5;

    /*
     * Phase 3: find e_sewn in [e_zero, e_extr].
     *
     * Sewing condition: d/dT[dE/dx] - dE/dx/T = 0.
     * Approximated numerically as a centred finite difference of dE/dx
     * with step h, minus dE/dx / T.  The golden-section search drives this
     * expression to zero.
     */
    e_min = e_zero;
    e_max = e_extr;
    rla = e_min + GOLD_LOWER * (e_max - e_min);
    rmu = e_min + GOLD_UPPER * (e_max - e_min);
    tla = (bethe_raw((rla + h) * a, proj, target) - bethe_raw((rla - h) * a, proj, target)) / h
          - bethe_raw(rla * a, proj, target) / rla;
    tmu = (bethe_raw((rmu + h) * a, proj, target) - bethe_raw((rmu - h) * a, proj, target)) / h
          - bethe_raw(rmu * a, proj, target) / rmu;

    while ((e_max - e_min) >= epsilon) {
        if (tla <= 0.0) {
            e_max = rmu;
            rmu = rla;
            tmu = tla;
            rla = e_min + GOLD_LOWER * (e_max - e_min);
            tla = (bethe_raw((rla + h) * a, proj, target) - bethe_raw((rla - h) * a, proj, target)) / h
                  - bethe_raw(rla * a, proj, target) / rla;
        } else {
            e_min = rla;
            rla = rmu;
            tla = tmu;
            rmu = e_min + GOLD_UPPER * (e_max - e_min);
            tmu = (bethe_raw((rmu + h) * a, proj, target) - bethe_raw((rmu - h) * a, proj, target)) / h
                  - bethe_raw(rmu * a, proj, target) / rmu;
        }
    }
    sewn->e_sewn = (e_min + e_max) * 0.5;
    sewn->f_sewn = bethe_raw(sewn->e_sewn * a, proj, target);
}

double osh_physics_bethe_eval(double t_per_nucleon,
                              struct osh_physics_bethe_projectile const *proj,
                              struct osh_physics_bethe_target const *target,
                              struct osh_physics_bethe_sewn const *sewn) {
    double t_total; /* total kinetic energy of the projectile [MeV] */

    t_total = t_per_nucleon * proj->a;

    /*
     * Lindhard-Scharff low-energy extension.
     *
     * Below the sewing point the Bethe formula loses validity (nuclear
     * stopping becomes dominant and the Born approximation breaks down).
     * The LSS formula dE/dx ~ sqrt(T) is used instead [LSS].  The constant
     * C_lss = f_sewn / sqrt(e_sewn * a_proj) ensures continuity at e_sewn.
     *
     * Note: the square-root dependence on *total* kinetic energy T = T/u * A,
     * so sqrt(T) = sqrt(e_per_nucleon * A), matching the libdedx convention.
     */
    if (t_per_nucleon <= sewn->e_sewn) {
        return (sewn->f_sewn / sqrt(sewn->e_sewn * proj->a)) * sqrt(t_total);
    }

    return bethe_raw(t_total, proj, target);
}

/* ---- Public: Hubert effective charge Z_eff -------------------------------- */

double osh_physics_bethe_z_eff(double t_per_nucleon, double proj_z, double proj_a, double target_z_mean) {
    double dd;
    double u1;
    double u2;
    double u3;
    double u4;
    double gamma_eff;

    if (target_z_mean <= 0.0 || proj_z <= 0.0) {
        return proj_z; /* guard: treat as fully stripped */
    }

    /*
     * Hubert effective-charge GAMMA [H89], identical to the calculation
     * inside bethe_raw().  t_per_nucleon corresponds to t_total_mev / proj->a
     * in that function.
     */
    dd = 1.164 + 0.2319 * exp(-0.004302 * target_z_mean);
    u1 = dd + 1.658 * exp(-0.05170 * proj_z);
    u2 = 8.144 + 0.09876 * log(target_z_mean);
    u3 = 0.3140 + 0.01072 * log(target_z_mean);
    u4 = 0.5218 + 0.02521 * log(target_z_mean);
    gamma_eff = 1.0 - u1 * exp((-u2 * pow(t_per_nucleon, u3)) / pow(proj_z, u4));

    if (gamma_eff < 0.0) {
        gamma_eff = 0.0;
    } else if (gamma_eff > 1.0) {
        gamma_eff = 1.0;
    }

    (void) proj_a; /* A not needed by the Hubert formula; kept for caller symmetry */
    return proj_z * gamma_eff;
}

/* ---- Internal Bethe evaluator -------------------------------------------- */

/**
 * @brief Evaluate the modified Bethe-Bloch formula at a given total kinetic energy.
 *
 * @details
 * Returns the mass stopping power [MeV cm²/g] at total kinetic energy
 * @p t_total_mev [MeV] (not per nucleon).  No low-energy clamping is applied;
 * the result can be negative below the Bragg-peak zero crossing.
 * Callers should use osh_physics_bethe_eval() which applies the LSS extension.
 *
 * The formula (in SI-adjacent mixed units matching [BB]):
 *
 *   -dE/dx = K * (Z_t/A_t) * (Z_eff^2 / beta^2)
 *            * [ ln(2 m_e c^2 beta^2 gamma^2 / I) - beta^2 - delta/2 ]
 *
 * where K = 0.307075 MeV cm^2/mol, I is the mean excitation energy [eV],
 * Z_eff = Z_proj * GAMMA (Hubert effective charge factor), and delta is the
 * Sternheimer density correction.
 *
 * @par References
 * [BB]  Bethe, Ann. Phys. 5, 325 (1930); Bloch, ibid. 16, 285 (1933). \n
 * [SP]  Sternheimer & Peierls, Phys. Rev. B 3, 3681 (1971). \n
 * [H89] Hubert, Bimbot & Gauvin, NIMB 36, 357 (1989).
 *
 * @param[in] t_total_mev  Total kinetic energy [MeV] (= T/u * A_proj).
 * @param[in] proj         Projectile parameters.
 * @param[in] target       Target parameters.
 */
static double bethe_raw(double t_total_mev,
                        struct osh_physics_bethe_projectile const *proj,
                        struct osh_physics_bethe_target const *target) {
    double mass;      /* projectile rest mass [MeV/c^2], taken from proj->mass_mev */
    double momentum;  /* projectile momentum [MeV/c] */
    double e0;        /* total energy [MeV] = T + m */
    double gamma2;    /* Lorentz factor squared: gamma^2 = (E/m)^2 */
    double beta2;     /* velocity squared: beta^2 = (gamma^2-1)/gamma^2 */
    double g2b2;      /* gamma^2 * beta^2 = (p/m)^2, used in Bethe log */
    double x;         /* Sternheimer log10(p/m) */
    double hnp;       /* plasma frequency of target [eV] */
    double c_stern;   /* Sternheimer material constant C */
    double cav;       /* cav = -C */
    double xal;       /* xal = cav / ln(10) */
    double x0;        /* lower Sternheimer threshold */
    double x1;        /* upper Sternheimer threshold */
    double alit;      /* cubic interpolation coefficient for Sternheimer */
    double delta;     /* density effect correction (dimensionless) */
    double dd;        /* Hubert Z_eff intermediate: target Z dependence */
    double u1;        /* Hubert parameter U1 */
    double u2;        /* Hubert parameter U2 */
    double u3;        /* Hubert parameter U3 */
    double u4;        /* Hubert parameter U4 */
    double gamma_eff; /* Hubert GAMMA: fractional effective charge */
    double z_eff2;    /* Z_eff^2 = Z_proj^2 * GAMMA^2 */
    double dedx;      /* result: mass stopping power [MeV cm^2/g] */

    mass = proj->mass_mev;

    /* ---- Relativistic kinematics ---------------------------------------- */

    /*
     * p = sqrt(T*(T + 2m)),  derived from E^2 = p^2 + m^2  with E = T + m.
     * Avoids catastrophic cancellation at low T << m.
     */
    momentum = sqrt(t_total_mev * (t_total_mev + 2.0 * mass));
    e0 = t_total_mev + mass;
    gamma2 = (e0 / mass) * (e0 / mass);
    beta2 = (gamma2 - 1.0) / gamma2;

    /*
     * gamma^2 * beta^2 = p^2 / m^2.  This is the kinematic factor that
     * appears in the argument of the Bethe logarithm:
     *   ln(2 m_e c^2 * gamma^2 beta^2 / I).
     */
    g2b2 = gamma2 * beta2;

    /* ---- Sternheimer-Peierls density effect correction [SP] -------------- */

    /*
     * The density effect correction delta reduces the stopping power at high
     * velocities due to the polarisation of the medium by the projectile's
     * electric field (relativistic rise suppression).
     *
     * x = log10(p/m) = log10(beta * gamma).
     */
    x = log10(momentum / mass);

    /*
     * Plasma frequency of the target medium [eV]:
     *   hnp = 28.8 * sqrt(rho * Z_target / A_target)
     *
     * 28.8 eV corresponds to 4 pi N_A r_e^(1/2) * (m_e c^2)^(1/2) in
     * appropriate units.  See [SP] eq. (3).
     */
    hnp = 28.8 * sqrt(target->rho * target->z_mean / target->a_mean);

    /*
     * Material constant C = -2 ln(I/hnp) - 1.
     * At high momenta delta → 2 ln(beta*gamma) + C, recovering the
     * classical density correction.
     */
    c_stern = -2.0 * log(target->i_value / hnp) - 1.0;
    cav = -c_stern;
    xal = cav / BETHE_LN10;

    /*
     * Sternheimer threshold parameters x0 and x1.
     * Two regimes depending on whether I >= 100 eV (conductors/heavy targets)
     * or I < 100 eV (gases and light materials such as hydrogen).
     * Values from [SP] Table II.
     */
    if (target->i_value >= 100.0) {
        x1 = 3.0;
        x0 = (cav >= 5.215) ? (0.326 * cav - 1.5) : 0.2;
    } else {
        x1 = 2.0;
        x0 = (cav >= 3.681) ? (0.326 * cav - 1.0) : 0.2;
    }

    /*
     * Cubic interpolation coefficient for the transition region x0 < x < x1.
     * Ensures continuity and a smooth first derivative of delta at x0 and x1.
     */
    alit = BETHE_LN10 * (xal - x0) / ((x1 - x0) * (x1 - x0) * (x1 - x0));

    /*
     * Piecewise definition of delta [SP] eq. (5):
     *   x <= x0              : delta = 0  (no correction at low momenta)
     *   x0 < x <= x1         : delta = ln(10)*x + c + alit*(x1-x)^3  (cubic)
     *   x > x1               : delta = ln(10)*x + c  (full relativistic rise)
     */
    if (x <= x0) {
        delta = 0.0;
    } else if (x <= x1) {
        delta = BETHE_LN10 * x + c_stern + alit * (x1 - x) * (x1 - x) * (x1 - x);
    } else {
        delta = BETHE_LN10 * x + c_stern;
    }

    /* ---- Hubert effective charge Z_eff [H89] ----------------------------- */

    /*
     * Heavy ions at clinical energies are not fully stripped; they capture
     * electrons from the target, reducing their effective charge below Z_proj.
     * The Hubert et al. parameterisation gives the fractional charge GAMMA
     * as a function of both projectile Z and target Z, fitted to measured data.
     *
     * GAMMA = 1 - U1 * exp(-U2 * (T/A)^U3 / Z_proj^U4)
     *
     * where T/A is kinetic energy per nucleon [MeV/nucleon] and the U_i are
     * target-Z-dependent coefficients from Table 1 of [H89].
     *
     * DD encodes the target-Z dependence of the high-energy limit of Z_eff.
     * U1 adds a projectile-Z correction.  U2, U3, U4 control the energy
     * dependence of the charge-stripping transition.
     *
     * For protons (Z_proj = 1) GAMMA ≈ 1 at all relevant energies, recovering
     * the standard Bethe formula.
     */
    dd = 1.164 + 0.2319 * exp(-0.004302 * target->z_mean);
    u1 = dd + 1.658 * exp(-0.05170 * proj->z);
    u2 = 8.144 + 0.09876 * log(target->z_mean);
    u3 = 0.3140 + 0.01072 * log(target->z_mean);
    u4 = 0.5218 + 0.02521 * log(target->z_mean);
    gamma_eff = 1.0 - u1 * exp((-u2 * pow(t_total_mev / proj->a, u3)) / pow(proj->z, u4));
    z_eff2 = proj->z * proj->z * gamma_eff * gamma_eff;

    /* ---- Modified Bethe-Bloch formula [BB] ------------------------------- */

    /*
     * -dE/dx = K * (Z_t/A_t) * (Z_eff^2/beta^2)
     *          * [ ln(2 m_e c^2 * beta^2 gamma^2 / I) - beta^2 - delta/2 ]
     *
     * The argument of the logarithm:
     *   (2 m_e c^2 / I) * gamma^2 beta^2
     *   = (1.022e6 eV / I[eV]) * g2b2
     *
     * Units: K [MeV cm^2/mol] * (Z/A) [mol/g] = [MeV cm^2/g].
     *
     * Shell corrections (C/Z term from [BB]) are not included here; both
     * libdedx and SHIELD-HIT omit them, relying on the ICRU-tabulated I-values
     * which partially absorb the shell correction effect.
     */
    dedx = BETHE_K * (target->z_mean / target->a_mean) * (z_eff2 / beta2)
           * (log((BETHE_2MEC2_EV / target->i_value) * g2b2) - beta2 - delta * 0.5);

    return dedx;
}
