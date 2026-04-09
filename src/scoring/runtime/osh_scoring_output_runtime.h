#ifndef OSH_SCORING_OUTPUT_RUNTIME_H
#define OSH_SCORING_OUTPUT_RUNTIME_H

#include <stddef.h>

#include "scoring/runtime/osh_scoring_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

struct osh_scoring_page_filter_ref {
    size_t filter_idx;
};

struct osh_scoring_page_settings_ref {
    size_t settings_idx;
};

struct osh_scoring_page_runtime {
    char *quantity;
    struct osh_scoring_page_filter_ref *filters;
    struct osh_scoring_page_settings_ref *settings;
    double *data;
    double *data_var;
    double *data2;
    double *data2_var;
    size_t output_idx;
    size_t geometry_idx;
    size_t nfilters;
    size_t nsettings;
    size_t len;
    enum osh_scoring_score_kind score_kind;
    char has_data2;
    char variance;
    char divide;
    char postproc;
};

struct osh_scoring_output_runtime {
    char *filename;
    char *fileformat;
    size_t geometry_idx;
    size_t *page_indices;
    size_t npages;
};

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_OUTPUT_RUNTIME_H */
