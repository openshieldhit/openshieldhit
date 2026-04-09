#ifndef OSH_SCORING_FILTER_RUNTIME_H
#define OSH_SCORING_FILTER_RUNTIME_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

struct osh_scoring_filter_runtime_rule {
    char field[16];
    char op[4];
    double value;
};

struct osh_scoring_filter_runtime {
    char *name;
    struct osh_scoring_filter_runtime_rule *rules;
    size_t nrules;
};

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_FILTER_RUNTIME_H */
