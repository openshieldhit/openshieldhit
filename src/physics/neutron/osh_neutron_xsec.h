#ifndef OSH_NEUTRON_XSEC_H
#define OSH_NEUTRON_XSEC_H

/**
 * @file osh_neutron_xsec.h
 * @brief Neutron cross-section lookup with JEFF-4.0 Tier-1 and optical Tier-2 fallback.
 *
 * @details
 * Tier-1: log-log interpolation on condensed JEFF-4.0 PENDF0K tables
 * (tools/condense_neutron_xsec.py).  The current grid has 31 energy points,
 * i.e. 30 interpolation intervals.  Channels: σ_tot, σ_el, σ(n,n'),
 * σ(n,2n), σ(n,γ), σ(n,p), σ(n,α).
 *
 * Tier-2 (any (Z,A) not in the JEFF tables): σ_R from Tripathi (1999),
 * σ_el from the geometric approximation π r₀² A^(2/3) (r₀ = 1.2 fm).
 * All sub-channels are zero; the reaction layer treats the non-elastic part
 * as a generic compound event.  A one-time diagnostic is emitted per missing
 * (Z,A) the first time it is encountered.
 *
 * The module is stateless beyond the one-time warning tracker; no dynamic
 * allocation is performed.  Call osh_neutron_xsec_compile() once at startup
 * and osh_neutron_xsec_free() at teardown.
 */

#include "openshieldhit/status.h"

#ifdef __cplusplus
extern "C" {
#endif

struct osh_diag_sink;

/** Maximum number of distinct Tier-2 (Z,A) fallbacks that will be logged. */
#define OSH_NEUTRON_XSEC_MAX_WARNED 64

/**
 * @brief Cross-section model state.
 *
 * Stack-allocatable.  The JEFF table is a static constant embedded in the
 * translation unit; this struct carries only the diagnostic sink and the
 * one-time fallback warning tracker.
 */
struct osh_neutron_xsec {
    struct osh_diag_sink const *diag;
    int warned_z[OSH_NEUTRON_XSEC_MAX_WARNED]; /**< (Z,A) pairs already warned */
    int warned_a[OSH_NEUTRON_XSEC_MAX_WARNED];
    int n_warned;
};

/**
 * @brief Cross-section values for one nuclide at one energy [mb].
 *
 * All seven JEFF channels plus the non-elastic sum are filled by
 * osh_neutron_xsec_lookup().
 */
struct osh_neutron_xsec_result {
    double tot;   /**< σ_total              (MT=1)   */
    double el;    /**< σ_elastic            (MT=2)   */
    double nn;    /**< σ(n,n') inelastic    (MT=4)   */
    double n2n;   /**< σ(n,2n)              (MT=16)  */
    double ng;    /**< σ(n,γ) capture       (MT=102) */
    double np;    /**< σ(n,p)               (MT=103) */
    double na;    /**< σ(n,α)               (MT=107) */
    double nonel; /**< σ_non-elastic = σ_tot − σ_el  */
};

/**
 * @brief Initialise the cross-section model.
 *
 * Stores the diagnostic sink and resets the fallback warning counter.  The
 * JEFF table is accessed as a static constant; no allocation is performed.
 *
 * @param[in]  diag  Diagnostic sink; may be NULL.
 * @param[out] xsec  Caller-allocated struct to initialise.
 * @returns OSH_OK, or OSH_EINVAL if xsec is NULL.
 */
enum osh_status osh_neutron_xsec_compile(struct osh_diag_sink const *diag, struct osh_neutron_xsec *xsec);

/**
 * @brief Release any resources held by the model.
 *
 * Currently a no-op (all storage is static), but pairs symmetrically with
 * compile() to follow the module lifecycle convention.
 */
void osh_neutron_xsec_free(struct osh_neutron_xsec *xsec);

/**
 * @brief Look up all cross-section channels for nuclide (Z, A) at energy E.
 *
 * Searches the JEFF-4.0 Tier-1 table first.  If the nuclide is absent, falls
 * back to Tripathi + geometric (Tier-2), emitting a one-time diagnostic.
 *
 * Binary search on the non-equidistant log-spaced energy grid; log-log
 * interpolation in each bracket.  Linear interpolation is used at threshold
 * crossings (where one endpoint is zero).  Clamped at the grid edges.
 *
 * @param[in,out] xsec   Model state (warning tracker is mutated on first fallback).
 * @param[in]     z      Target atomic number.
 * @param[in]     a      Target mass number.
 * @param[in]     e_mev  Neutron kinetic energy [MeV] (lab frame, target at rest).
 * @param[out]    out    Filled with σ values [mb] and σ_nonel.
 */
void osh_neutron_xsec_lookup(
    struct osh_neutron_xsec *xsec, int z, int a, double e_mev, struct osh_neutron_xsec_result *out);

#ifdef __cplusplus
}
#endif

#endif /* OSH_NEUTRON_XSEC_H */
