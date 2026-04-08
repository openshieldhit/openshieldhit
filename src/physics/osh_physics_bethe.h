#ifndef OSH_PHYSICS_BETHE_H
#define OSH_PHYSICS_BETHE_H

/**
 * @file osh_physics_bethe.h
 *
 * @brief Modified Bethe-Bloch mass stopping power for ions in matter.
 *
 * @details
 * Implements the modified Bethe-Bloch formula with:
 *  - Sternheimer-Peierls density effect correction (delta)
 *  - Hubert et al. effective-charge parameterisation (Z_eff)
 *  - Lindhard-Scharff low-energy extension, stitched at the sewing point
 *
 * The formulas, coefficients, and low-energy sewing procedure are taken
 * directly from the primary literature listed below.  Numerical values were
 * cross-checked for consistency against the libdedx library (Herrmann et al.,
 * https://github.com/nbassler/libdedx, GPL-3.0), which implements the same
 * physics for the same application domain.
 *
 * @par References
 * [BB]  Bethe, Ann. Phys. 5, 325 (1930); Bloch, ibid. 16, 285 (1933). \n
 * [SP]  Sternheimer & Peierls, Phys. Rev. B 3, 3681 (1971). Density effect. \n
 * [H89] Hubert, Bimbot & Gauvin, NIMB 36, 357 (1989). Effective charge. \n
 * [LSS] Lindhard, Scharff & Schiott, Mat. Fys. Medd. 33:14 (1963). \n
 *       Low-energy nuclear+electronic stopping.
 */

/**
 * @brief Target material parameters needed by the Bethe model.
 *
 * @details
 * Separate z_mean and a_mean fields match the libdedx convention (TZ0, TA0)
 * (which matches SHIELD-HIT).  In the formula they always appear as the ratio
 * Z/A; a future refactor could replace both with a single z_over_a field:
 *
 *   z_over_a = sum_i(w_i * Z_i / A_i)   [mol/g, mass-fraction weighted]
 *
 * That would make the compound case unambiguous and remove the temptation to
 * pass per-element Z and A independently.  For now the two-field layout is
 * kept for direct traceability to both reference implementations.
 *
 * For a compound, the caller must supply effective values weighted by mass
 * fraction: z_mean = sum_i(w_i * Z_i),  a_mean such that z_mean/a_mean
 * = sum_i(w_i * Z_i / A_i).
 */
struct osh_physics_bethe_target {
    double z_mean;  /* effective atomic number (mass-fraction weighted) */
    double a_mean;  /* effective atomic mass such that z_mean/a_mean = sum_i(w_i*Z_i/A_i) [Da] */
    double rho;     /* density [g/cm³] */
    double i_value; /* mean excitation energy [eV] */
};

/**
 * @brief Projectile parameters needed by the Bethe model.
 *
 * @details
 * The kinematics use the fully-stripped nuclear rest mass supplied in
 * @p mass_mev. Callers should obtain it from the particle module, e.g.
 * via osh_particle_nuclear_mass_mev_from_za(), rather than using the old
 * 940*A approximation.
 */
struct osh_physics_bethe_projectile {
    double z;        /* atomic number */
    double a;        /* integer nucleon number A */
    double mass_mev; /* fully-stripped nuclear rest mass [MeV/c²] */
};

/**
 * @brief Pre-computed Lindhard-Scharff sewing point for one (proj, target) pair.
 *
 * @details
 * The Bethe formula gives unphysical (negative or rising) results below a
 * material- and projectile-dependent threshold.  The Lindhard-Scharff formula
 *   dE/dx = (f_sewn / sqrt(e_sewn * a_proj)) * sqrt(T)
 * is used instead below the sewing energy.  The sewing point is the kinetic
 * energy per nucleon [MeV/nucleon] where the two curves are tangentially joined
 * (slope continuity), found once by golden-section search at setup time.
 */
struct osh_physics_bethe_sewn {
    double e_sewn; /* sewing energy per nucleon [MeV/nucleon] */
    double f_sewn; /* Bethe dE/dx at sewing point [MeV cm²/g] */
};

/**
 * @brief Find the Lindhard-Scharff sewing point for a given (proj, target) pair.
 *
 * @details
 * This is a one-time setup call per (projectile, material) combination.  The
 * result should be cached and passed to every osh_physics_bethe_eval() call for
 * the same pair.  The golden-section search is identical to the algorithm in
 * libdedx gold_section().
 *
 * @par References
 * [LSS] Lindhard, Scharff & Schiott, Mat. Fys. Medd. 33:14 (1963). \n
 *       Low-energy extension. \n
 * libdedx `gold_section()` implementation for the numerical sewing procedure.
 *
 * @param[in]  proj    Projectile parameters.
 * @param[in]  target  Target material parameters.
 * @param[out] sewn    Receives the sewing energy and stopping power.
 */
void osh_physics_bethe_sewn_compute(struct osh_physics_bethe_projectile const *proj,
                                    struct osh_physics_bethe_target const *target,
                                    struct osh_physics_bethe_sewn *sewn);

/**
 * @brief Evaluate mass stopping power at a single kinetic energy.
 *
 * @details
 * Returns the Bethe-Bloch result above the sewing point, and the
 * Lindhard-Scharff sqrt(T) extension below it.  The @p sewn struct must have
 * been filled by osh_physics_bethe_sewn_compute() for this (proj, target) pair.
 *
 * @par References
 * [BB]  Bethe, Ann. Phys. 5, 325 (1930); Bloch, ibid. 16, 285 (1933). \n
 * [SP]  Sternheimer & Peierls, Phys. Rev. B 3, 3681 (1971). \n
 * [H89] Hubert, Bimbot & Gauvin, NIMB 36, 357 (1989). \n
 * [LSS] Lindhard, Scharff & Schiott, Mat. Fys. Medd. 33:14 (1963).
 *
 * @param[in] t_per_nucleon  Kinetic energy per nucleon [MeV/nucleon].
 * @param[in] proj           Projectile parameters.
 * @param[in] target         Target material parameters.
 * @param[in] sewn           Pre-computed sewing point.
 *
 * @returns Mass stopping power [MeV cm²/g], always >= 0.
 */
double osh_physics_bethe_eval(double t_per_nucleon,
                              struct osh_physics_bethe_projectile const *proj,
                              struct osh_physics_bethe_target const *target,
                              struct osh_physics_bethe_sewn const *sewn);

#endif /* OSH_PHYSICS_BETHE_H */
