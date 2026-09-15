#ifndef OSH_SCORING_MCPL_RECORD_H
#define OSH_SCORING_MCPL_RECORD_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief One phase-space particle booked by the MCPL estimator (issue #328).
 *
 * @details
 * Filled by @c osh_scoring_estimator_step_mcpl() in openshieldhit's own units —
 * position in cm, energy in MeV — which are already MCPL's units, so
 * @c osh_scoring_save_mcpl_output() converts none of them; it only packs this
 * struct into an @c mcpl_particle_t.  Kept independent of the vendored
 * <mcpl.h> so the hot-path headers never pull in a third-party dependency;
 * only the save layer includes <mcpl.h>.
 */
struct osh_scoring_mcpl_record {
    double position[3];  /* [cm] */
    double direction[3]; /* unit vector */
    double ekin;         /* kinetic energy [MeV] */
    double weight;       /* st->wt at the moment this record was booked */
    int32_t pdgcode;     /* PDG Monte Carlo particle numbering scheme code */
    uint32_t gen;        /* st->gen: 0 = beam primary, >=1 = secondary generation */
};

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_MCPL_RECORD_H */
