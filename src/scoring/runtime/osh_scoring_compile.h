#ifndef OSH_SCORING_COMPILE_H
#define OSH_SCORING_COMPILE_H

#include "openshieldhit/diag.h"
#include "openshieldhit/status.h"
#include "scoring/osh_scoring.h"
#include "scoring/runtime/osh_scoring_runtime.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Compile parsed scorer configuration into scorer runtime objects.
 *
 * @details
 * This resolves raw `detect.dat` names to dense indices, computes geometry bin
 * counts, groups pages by shared geometry, and allocates page-local scoring
 * buffers. Geometry traversal itself is not built here yet; this is the first
 * compiled bridge from readable parser structs to simulation-ready scorer
 * runtime objects.
 */
enum osh_status osh_scoring_compile(struct osh_scoring_workspace const *ws,
                                    struct osh_diag_sink const *diag,
                                    struct osh_scoring_runtime *rt);

/**
 * @brief Merge rt->settings[] into each page's embedded sset.
 *
 * @details
 * Must be called after any post-compile mutation of rt->settings (e.g. after
 * material-name resolution in the simulation layer sets has_medium).  Safe to
 * call multiple times; each call rebuilds sset from scratch.
 */
void osh_scoring_runtime_finalize_ssets(struct osh_scoring_runtime *rt);

/**
 * @brief Release a compiled scorer runtime and all owned memory.
 */
void osh_scoring_runtime_free(struct osh_scoring_runtime *rt);

/**
 * @brief Shape a caller-provided accumulator set to match @p rt's pages.
 *
 * @details
 * Fills the caller-owned, @c rt->npages-long array @p set so each element is sized
 * and shaped like the matching master page (same @c len, same @c data2 presence),
 * all zero-initialised.  This is the per-worker deposit target a replica (issue
 * #230) or a future parallel worker (#161/#195) scores into, then folds into the
 * master with @ref osh_scoring_accumulator_merge.  Because a page's accumulator and
 * its clone share the page descriptor, the merge always agrees on optional-array
 * presence.
 *
 * The caller owns @p set's storage (e.g. one slice of a flat per-replica block), so
 * this allocates only the per-page arrays — never the set array itself — and all
 * allocation happens at run setup, never on the hot path (DEVELOPER.md §10).
 *
 * @param[in]  rt   Compiled scoring runtime to mirror.
 * @param[out] set  Caller-owned array of at least @c rt->npages accumulators to
 *                  populate; released with @ref osh_scoring_runtime_free_accumulator_set.
 * @returns OSH_OK on success (including @c rt->npages == 0, a no-op), OSH_EINVAL on
 *          NULL argument, OSH_ENOMEM on allocation failure (the pages populated so
 *          far are freed, so @p set is left all-zeroed and leaks nothing).
 */
enum osh_status osh_scoring_runtime_alloc_accumulator_set(struct osh_scoring_runtime const *rt,
                                                          struct osh_scoring_accumulator *set);

/**
 * @brief Free the per-page arrays of a set from @ref osh_scoring_runtime_alloc_accumulator_set.
 *
 * @details
 * Frees each accumulator's arrays and zeroes it, but NOT the @p set array itself
 * (the caller owns that storage).  Safe with a NULL @p set (no-op).  @p npages must
 * be the count the set was shaped with.
 */
void osh_scoring_runtime_free_accumulator_set(struct osh_scoring_accumulator *set, size_t npages);

/**
 * @brief Allocate a private traversal scratch sized like @p rt's master scratch.
 *
 * @details
 * Gives a worker its own @ref osh_scoring_scratch with a @c crossing_buf pre-sized
 * to the largest geometry (mirroring @c rt->master_scratch), so the deposit path
 * never grows it on the hot path — @ref osh_scoring_score_step rejects an
 * under-sized scratch rather than reallocating.  A runtime with no crossing
 * geometry yields a NULL buffer (cap 0), exactly like the master.
 *
 * @param[in]  rt           Compiled scoring runtime to mirror.
 * @param[out] scratch_out  Receives the scratch (owned by the caller; free with
 *                          @ref osh_scoring_runtime_free_scratch).
 * @returns OSH_OK on success, OSH_EINVAL on NULL argument, OSH_ENOMEM on failure.
 */
enum osh_status osh_scoring_runtime_clone_scratch(struct osh_scoring_runtime const *rt,
                                                  struct osh_scoring_scratch *scratch_out);

/**
 * @brief Free a private scratch from @ref osh_scoring_runtime_clone_scratch.
 *
 * Releases @c crossing_buf and zeroes the struct.  Safe with a NULL @p scratch.
 */
void osh_scoring_runtime_free_scratch(struct osh_scoring_scratch *scratch);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_COMPILE_H */
