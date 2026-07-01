#ifndef OSH_PHYSICS_SCAT_H
#define OSH_PHYSICS_SCAT_H

struct osh_rng;

/** MCS angular model; integer values mirror enum osh_transport_mcs_mode. */
enum osh_mcs_model {
    OSH_MCS_OFF = 0,
    OSH_MCS_GAUSSIAN = 1, /**< Highland Gaussian core (no tail). */
    OSH_MCS_MOLIERE = 2,  /**< Full Bethe-Moliere (core + Rutherford tail). */
    OSH_MCS_WENTZEL = 3   /**< Reserved; not implemented here. */
};

/**
 * @brief Sample and apply an MCS deflection for the selected model.
 *
 * @details
 * Dispatches on @p model.  The dispatcher owns no model-specific state:
 * Highland and Bethe-Moliere details live in their model modules, and future
 * Wentzel support should be added as another sibling module.
 *
 * @param[in]  model            MCS model.
 * @param[in]  v                Incident unit direction (length 3).
 * @param[out] w                Exit unit direction (length 3); always set.
 * @param[in]  t_kin_mev        Kinetic energy at mid-step [MeV].
 * @param[in]  mass_mev         Projectile rest mass [MeV/c^2].
 * @param[in]  z_eff            Effective projectile charge.
 * @param[in]  d_gcm2           Substep areal density rho*ds [g/cm^2].
 * @param[in]  path_scale_gcm2  Macroscopic path scale [g/cm^2].
 * @param[in]  x0_gcm2          Radiation length [g/cm^2] for Highland.
 * @param[in]  chic2_coeff      Medium chi_c^2 coefficient for Moliere.
 * @param[in]  screen_z         Medium effective screening Z for Moliere.
 * @param[in]  rng              RNG state.
 * @returns 1 if a deflection was applied, 0 if @p w was copied unchanged.
 */
int osh_physics_mcs_scatter(enum osh_mcs_model model,
                            double const v[3],
                            double w[3],
                            double t_kin_mev,
                            double mass_mev,
                            double z_eff,
                            double d_gcm2,
                            double path_scale_gcm2,
                            double x0_gcm2,
                            double chic2_coeff,
                            double screen_z,
                            struct osh_rng *rng);

#endif /* OSH_PHYSICS_SCAT_H */
