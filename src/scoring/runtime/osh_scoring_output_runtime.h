#ifndef OSH_SCORING_OUTPUT_RUNTIME_H
#define OSH_SCORING_OUTPUT_RUNTIME_H

#include <stddef.h>

#include "scoring/runtime/osh_scoring_defs.h"
#include "scoring/runtime/osh_scoring_filter_runtime.h"

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Index reference from a page to a compiled settings block. */
struct osh_scoring_page_settings_ref {
    size_t settings_idx; /* Index into osh_scoring_runtime::settings[]. */
};

/**
 * @brief Hot-path override fields condensed from a page's Settings block(s).
 *
 * Contains only the fields read by the scoring hot path.  Much smaller than
 * the full osh_scoring_settings_runtime, reducing page-struct cache pressure.
 */
struct osh_scoring_page_override {
    double density_g_cm3; /* Local density override [g/cm³]. */
    int medium;           /* Transport medium override index. */
    char has_medium;
    char has_density_g_cm3;
};

/**
 * @brief One compiled scoring page (a single quantity accumulator).
 *
 * @details
 * Each page owns a flat @p data array of length @p len = product of geometry
 * bin counts.  @p data2 / @p data2_var hold the secondary accumulator used by
 * two-pass LET averages; both are NULL for simple scorers.
 *
 * Index layout follows the geometry axis order, innermost-axis fastest
 * (row-major):  idx = ix + nx * (iy + ny * iz)
 */
struct osh_scoring_page_runtime {
    char *quantity;                                     /* Quantity keyword, lowercase (owned). */
    struct osh_scoring_filter_runtime_rule *flat_rules; /* Flattened filter rules, all ANDed (owned). */
    struct osh_scoring_page_settings_ref *settings;     /* Settings index references (owned). */
    double *data;                                       /* Primary accumulator array (owned). */
    double *data_var;                                   /* Variance accumulator for data (owned, may be NULL). */
    double *data2;                          /* Secondary accumulator for LET averages (owned, may be NULL). */
    double *data2_var;                      /* Variance for data2 (owned, may be NULL). */
    size_t output_idx;                      /* Index of the owning output in the runtime. */
    size_t geometry_idx;                    /* Index of the owning geometry in the runtime. */
    size_t nflat_rules;                     /* Number of flat filter rules. */
    size_t nsettings;                       /* Number of settings references. */
    size_t len;                             /* Total number of bins (= product of axis nbins). */
    enum osh_scoring_score_kind score_kind; /* What physical quantity is accumulated. */
    enum osh_scoring_postproc postproc;     /* How to combine across simulation runs. */
    struct osh_scoring_page_override sset;  /* Hot-path override fields, merged from all settings references. */
    char has_sset;                          /* Non-zero when at least one settings block is referenced. */
    char has_data2;                         /* Non-zero when data2/data2_var are allocated. */
    char variance;                          /* Non-zero when variance tracking is active. */
    char divide;                            /* Non-zero when bin values should be divided by bin volume. */
    /* Differential axis — all zero when no differential scoring (diff_nbins == 0). */
    size_t diff_nbins;                    /* Number of differential bins; 0 = plain scorer. */
    size_t diff_stride;                   /* Spatial bin count (= geo_nbins); stride per diff bin in data[]. */
    double diff_lo;                       /* Lower bound of the differential axis. */
    double diff_hi;                       /* Upper bound of the differential axis. */
    int diff_log;                         /* 0 = linear binning, 1 = logarithmic binning. */
    enum osh_scoring_diff_kind diff_kind; /* Physical quantity that determines the diff bin. */
    /* Second differential axis — all zero when unused (diff2_nbins == 0). */
    size_t diff2_nbins;  /* 0 = no second diff axis. */
    size_t diff2_stride; /* = diff_nbins * diff_stride; stride per diff2 bin in data[]. */
    double diff2_lo;
    double diff2_hi;
    int diff2_log;
    enum osh_scoring_diff_kind diff2_kind;
    /* Optional per-axis stopping-power override — e.g. "Diff1Type DEDX in_Si".
     * When has_diff_sset/has_diff2_sset is set and the Settings resolves to a
     * medium, LET and QEFF axis values are computed in that medium rather than
     * the transport medium.  Indices are into osh_scoring_runtime::settings[].
     * The resolved medium/density are cached in diff_sset/diff2_sset after
     * osh_scoring_runtime_finalize_ssets() runs for the second time. */
    size_t diff_sset_idx;                        /* index into rt->settings[]; valid when has_diff_sset */
    size_t diff2_sset_idx;                       /* index into rt->settings[]; valid when has_diff2_sset */
    struct osh_scoring_page_override diff_sset;  /* resolved override for diff1 axis */
    struct osh_scoring_page_override diff2_sset; /* resolved override for diff2 axis */
    char has_diff_sset;                          /* non-zero when diff1 axis has a settings override */
    char has_diff2_sset;                         /* non-zero when diff2 axis has a settings override */
};

/**
 * @brief Compiled output: one file with an ordered list of page indices.
 */
struct osh_scoring_output_runtime {
    char *filename;       /* Output file name (owned). */
    char *fileformat;     /* Format keyword, lowercase (owned; NULL = default BDO). */
    size_t geometry_idx;  /* Index of the associated geometry. */
    size_t *page_indices; /* Ordered page indices written to this file (owned). */
    size_t npages;        /* Number of pages. */
};

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_OUTPUT_RUNTIME_H */
