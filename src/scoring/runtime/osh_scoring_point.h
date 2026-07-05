#ifndef OSH_SCORING_POINT_H
#define OSH_SCORING_POINT_H

#include "common/osh_step.h"
#include "openshieldhit/status.h"
#include "particle/osh_particle.h"
#include "scoring/runtime/osh_scoring_accumulator.h"
#include "scoring/runtime/osh_scoring_runtime.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Score one point energy deposit into the compiled scoring runtime.
 *
 * @details
 * The counterpart to @ref osh_scoring_score_step for particles that deposit
 * their energy at a single point with no track length — a recoil or fragment
 * born below the transport threshold, or a neutron-reaction local deposit
 * (c.f. issue #179).
 *
 * The signature mirrors @ref osh_scoring_score_step and reuses the same
 * @ref step carrier: @c st->de is the energy to deposit [MeV] at the point
 * @c st->p, with @c st->rho the local density, @c st->medium / @c st->zone the
 * location, and @c st->gen / @c st->wt / @c st->p[3] the history attributes used
 * by page filters and differential axes.  @c st->q and @c st->ds are ignored
 * (a point has no exit point or track length); callers may set @c st->q == st->p.
 *
 * Only ENERGY and DOSE/DOSEGY pages receive a point contribution for now.
 * FLUENCE and LET/QEFF have no meaning without a track length and are wired in
 * later via precomputed LETd/LETt-versus-Ekin tables; those page kinds are
 * silently skipped. Mesh (X,Y,Z), Cyl (R,Z), and Zone geometries are supported.
 *
 * @param rt       Read-only compiled scoring descriptor.
 * @param acc_set  Mutable accumulator storage (indexed in lockstep with rt->pages).
 * @param scratch  Caller-owned scratch (unused today; kept for score_step parity).
 * @param part     Depositing particle species (for filter/attribution).
 * @param st       Point carrier: st->p the location, st->de the energy [MeV].
 */
enum osh_status osh_scoring_score_point(struct osh_scoring_runtime const *rt,
                                        struct osh_scoring_accumulator *acc_set,
                                        struct osh_scoring_scratch *scratch,
                                        struct particle const *part,
                                        struct step const *st);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_POINT_H */
