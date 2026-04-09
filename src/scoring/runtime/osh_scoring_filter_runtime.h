#ifndef OSH_SCORING_FILTER_RUNTIME_H
#define OSH_SCORING_FILTER_RUNTIME_H

#include <stddef.h>

#include "particle/osh_particle.h"
#include "scoring/runtime/osh_scoring_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

struct osh_scoring_page_runtime;
struct osh_scoring_runtime;

struct osh_scoring_filter_runtime_rule {
    enum osh_scoring_filter_field field;
    enum osh_scoring_filter_op op;
    double value;
};

struct osh_scoring_filter_runtime {
    char *name;
    struct osh_scoring_filter_runtime_rule *rules;
    size_t nrules;
};

/**
 * @brief Evaluate all compiled filters attached to one scoring page.
 *
 * @details
 * All attached filters and all rules within each filter are combined with
 * logical AND. The page passes only if every rule evaluates true.
 */
int osh_scoring_page_passes_filters(struct osh_scoring_runtime const *rt,
                                    struct osh_scoring_page_runtime const *page,
                                    struct particle const *part,
                                    double energy_mev);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_FILTER_RUNTIME_H */
