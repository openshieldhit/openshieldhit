#ifndef OSH_NEUTRON_CONST_H
#define OSH_NEUTRON_CONST_H

/**
 * @file osh_neutron_const.h
 * @brief Neutron-wide physics constants shared by the transport driver and the
 *        reaction/cross-section layer.
 *
 * @details
 * Single home for the handful of energy constants that both the transport loop
 * (osh_transport_neutron.c) and the reaction sampler (osh_neutron_reaction.c)
 * need, so no value is defined twice.
 *
 * Ordering invariant (do not break):
 *
 *   table floor (1 meV)  ≤  OSH_NEUTRON_CUTOFF_DEFAULT_MEV  <  OSH_NEUTRON_THERMAL_E_MEV
 *
 * The cross-section tables extend down to 1 meV (1e-9 MeV); the default
 * transport cutoff sits at that floor, and the thermal energy sits strictly
 * above it.  A neutron that has thermalised is pinned to OSH_NEUTRON_THERMAL_E_MEV
 * (see the freeze branch in do_elastic()), which is therefore always above the
 * cutoff — so a thermal neutron never trips the energy-cutoff kill and only
 * terminates on capture or geometry escape.
 */

/**
 * Room-temperature thermal neutron energy [MeV]: kT at 20 °C = 0.0253 eV.
 *
 * Below this energy elastic collisions change direction only (energy frozen):
 * a one-group thermal model that avoids the free-gas 0 K cross sections
 * down-scattering a neutron unphysically toward zero, and sidesteps S(α,β)
 * molecular up-scattering.  Neutrons crossing this threshold from above are
 * clamped up to it, so all thermal neutrons share this single energy.
 */
#define OSH_NEUTRON_THERMAL_E_MEV 2.53e-8

/**
 * Default lower neutron transport cutoff [MeV], used when NEUTRLCUT <= 0.
 *
 * Set to the cross-section table floor (1 meV) so neutrons transport all the
 * way to thermal by default.  NEUTRLCUT > 0 raises this (e.g. 1e-3 restores the
 * historical 1 keV kill, disabling thermal transport).  Must stay strictly
 * below OSH_NEUTRON_THERMAL_E_MEV (see the invariant above).
 */
#define OSH_NEUTRON_CUTOFF_DEFAULT_MEV 1.0e-9f

#endif /* OSH_NEUTRON_CONST_H */
