#ifndef OSH_SCORING_SHADOW_H
#define OSH_SCORING_SHADOW_H

#include <stddef.h>
#include <stdint.h>

#include "openshieldhit/status.h"
#include "scoring/runtime/osh_scoring_runtime.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief A non-destructive presentation view of a live scoring runtime.
 *
 * @details
 * Holds a @ref osh_scoring_runtime (@c view) whose page array is a struct-copy of
 * the live pages — so every array and scalar is *aliased* — except the
 * @c acc.data of pages whose postprocess writes data (DOSEGY, LET/Qeff), which is
 * redirected to a private scratch buffer.  An out-of-place postprocess then
 * writes only those scratch buffers, leaving the live accumulators byte-identical.
 *
 * Scratch is allocated **once** (lazily, on the first @ref osh_scoring_shadow_refresh)
 * and reused for every subsequent snapshot; it is freed only at
 * @ref osh_scoring_shadow_free.  This is what keeps periodic dumps free of per-dump
 * allocation churn (see issue #170 / #191).  The shadow owns only its @c pages
 * array and the scratch buffers — never the aliased live arrays.
 */
struct osh_scoring_shadow {
    struct osh_scoring_runtime view;        /* handed to postprocess_into + the sink */
    struct osh_scoring_runtime const *live; /* source runtime; never mutated */
    struct osh_scoring_page_runtime *pages; /* owned: npages struct-copies of the live pages */
    double **scratch;                       /* owned: per-page private data buffer (NULL = aliased) */
    size_t npages;
    int scratch_ready; /* non-zero once the scratch buffers have been allocated */
};

/**
 * @brief Bind a shadow to a live runtime (no scratch allocated yet).
 *
 * @returns OSH_OK, OSH_EINVAL on NULL args, OSH_ENOMEM on allocation failure.
 */
enum osh_status osh_scoring_shadow_init(struct osh_scoring_shadow *shadow, struct osh_scoring_runtime const *live);

/**
 * @brief Lazily allocate scratch (first call), then postprocess the live runtime
 *        into the view.  Never mutates the live accumulators.
 */
enum osh_status osh_scoring_shadow_refresh(struct osh_scoring_shadow *shadow);

/**
 * @brief Bytes of private scratch the shadow needs: sum of @c data sizes over the
 *        pages whose postprocess writes data.  Independent of whether the scratch
 *        is allocated yet.
 */
uint64_t osh_scoring_shadow_bytes(struct osh_scoring_shadow const *shadow);

/** @brief The presentation runtime to hand a sink after a refresh (NULL-safe). */
struct osh_scoring_runtime const *osh_scoring_shadow_view(struct osh_scoring_shadow const *shadow);

/** @brief Free the owned pages array + scratch buffers; never the aliased live arrays. */
void osh_scoring_shadow_free(struct osh_scoring_shadow *shadow);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_SHADOW_H */
