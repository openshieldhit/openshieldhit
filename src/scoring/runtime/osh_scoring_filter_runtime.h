#ifndef OSH_SCORING_FILTER_RUNTIME_H
#define OSH_SCORING_FILTER_RUNTIME_H

#include <stddef.h>

#include "common/osh_step.h"
#include "particle/osh_particle.h"
#include "scoring/runtime/osh_scoring_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

struct osh_scoring_page_runtime;
struct osh_scoring_runtime;

/**
 * @brief One compiled filter rule, ready for hot-path evaluation.
 *
 * @details
 * Resolved from the raw @ref osh_scoring_filter_rule string form during
 * @ref osh_scoring_compile.  Integer codes replace string comparisons so
 * the hot path can evaluate rules with a simple switch/if chain.
 */
struct osh_scoring_filter_runtime_rule {
    enum osh_scoring_filter_field field; /* Which particle property to test. */
    enum osh_scoring_filter_op op;       /* Comparison operator. */
    double value;                        /* Right-hand side of the comparison. */
};

/**
 * @brief Compiled filter: a named set of rules evaluated in AND combination.
 */
struct osh_scoring_filter_runtime {
    char *name;                                    /* Filter name (owned). */
    struct osh_scoring_filter_runtime_rule *rules; /* Rule array (owned). */
    size_t nrules;                                 /* Number of rules. */
};

/**
 * @brief Evaluate all compiled filters attached to one scoring page.
 *
 * @details
 * All attached filters and all rules within each filter are combined with
 * logical AND. The page passes only if every rule evaluates true.
 *
 * Species properties (Z, A, mass, PDG) are read from @p part.
 * Per-history properties (energy, generation, primary index, weight) are read
 * from @p st, which carries them as set by the transport engine.
 */
int osh_scoring_page_passes_filters(struct osh_scoring_runtime const *rt,
                                    struct osh_scoring_page_runtime const *page,
                                    struct particle const *part,
                                    struct step const *st);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_FILTER_RUNTIME_H */
