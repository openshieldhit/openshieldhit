#ifndef OSH_SCORING_SAVE_MCPL_H
#define OSH_SCORING_SAVE_MCPL_H

#include "openshieldhit/status.h"
#include "scoring/save/osh_scoring_save.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Save one compiled output as an MCPL phase-space file (issue #328).
 *
 * @details
 * Serialises the page's mcpl_records[0..*mcpl_count) — booked by
 * osh_scoring_estimator_step_mcpl() during transport — into a real .mcpl file
 * via the vendored MCPL writer (src/thirdparty/mcpl). This is the only point
 * anywhere in the MCPL path that calls into the MCPL library: score_step only
 * ever touches openshieldhit's own osh_scoring_mcpl_record buffer.
 *
 * The file stores double-precision fields, per-particle userflags carrying the
 * generation number (0 = beam primary; documented in a header comment so a
 * reader knows the convention), and a "nstat" stat:sum: header entry recording
 * @p nstat — the total number of primaries this run represents, letting a
 * downstream consumer scale a partial phase-space dump back to a per-primary
 * basis the way a restart run needs to.
 *
 * Only single-page outputs are supported (osh_scoring_compile() already
 * rejects any other shape for FileFormat MCPL); this returns OSH_ENOTSUP if
 * that invariant is somehow violated, and OSH_ESTATE if the page is not an
 * MCPL page.
 *
 * @param[in] ws          Scoring workspace with output metadata and file paths.
 * @param[in] rt          Compiled scoring runtime with the booked records.
 * @param[in] nstat       Actual number of primary particles simulated; must be > 0.
 * @param[in] output_idx  Zero-based index into ws->outputs / rt->outputs.
 */
enum osh_status osh_scoring_save_mcpl_output(struct osh_scoring_workspace const *ws,
                                             struct osh_scoring_runtime const *rt,
                                             unsigned long long nstat,
                                             size_t output_idx);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_SAVE_MCPL_H */
