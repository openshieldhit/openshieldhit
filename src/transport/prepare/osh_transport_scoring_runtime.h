#ifndef OSH_TRANSPORT_SCORING_RUNTIME_H
#define OSH_TRANSPORT_SCORING_RUNTIME_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief One compiled scoring filter reference attached to a runtime page.
 *
 * @details
 * The raw parser keeps filter names. The later scoring-prepare step can
 * resolve those names into dense indices and precompiled predicates.
 */
struct osh_transport_scoring_page_filter_ref {
    size_t filter_idx; /* Dense runtime filter index. */
};

/**
 * @brief One compiled scoring settings reference attached to a runtime page.
 *
 * @details
 * Like filters, settings are parsed by name and later resolved into dense
 * runtime indices during the scoring-prepare/finalize step.
 */
struct osh_transport_scoring_page_settings_ref {
    size_t settings_idx; /* Dense runtime settings index. */
};

/**
 * @brief One compiled scoring page/quantity attached to a shared geometry.
 *
 * @details
 * This is the transport-facing form of one BDO page. Multiple pages may share
 * one geometry runtime object so Siddon/Jacobs traversal runs only once per
 * scored step before all page accumulators are updated.
 */
struct osh_transport_scoring_page_runtime {
    char *quantity; /* Quantity name, e.g. "DOSE" or "FLUENCE". */
    char *label;    /* Optional user-facing page label or output name. */

    struct osh_transport_scoring_page_filter_ref *filters;    /* Resolved filter refs. */
    struct osh_transport_scoring_page_settings_ref *settings; /* Resolved settings refs. */
    size_t nfilters;                                          /* Length of @ref filters. */
    size_t nsettings;                                         /* Length of @ref settings. */

    double *data;      /* Primary page accumulator. */
    double *data_var;  /* Optional variance for @ref data. */
    double *data2;     /* Optional second page accumulator. */
    double *data2_var; /* Optional variance for @ref data2. */
    size_t len;        /* Length of all active data arrays. */

    unsigned int quantity_id; /* Future dense quantity enum. */
    unsigned char has_data2;  /* Whether @ref data2 / @ref data2_var are allocated. */
    unsigned char variance;   /* Whether variance arrays are allocated. */
    unsigned char divide;     /* Whether post-save result should divide data by data2. */
    unsigned char postproc;   /* Future postprocessing mode enum/flag. */
};

/**
 * @brief One compiled scoring geometry with many attached pages.
 *
 * @details
 * The shared geometry object is the important performance invariant carried
 * over from SHIELD-HIT: one geometry traversal, many page updates.
 *
 * Future prepare/finalize work will resolve the raw geometry definition into:
 * - a geometry kind enum
 * - precomputed bin strides
 * - inverse volumes
 * - Siddon/Jacobs workspace
 */
struct osh_transport_scoring_geometry_runtime {
    char *name;  /* Geometry name. */
    char *kind;  /* Geometry kind, e.g. "Mesh", "Cyl", "Zone". */
    char *fname; /* Output file name, typically one BDO file per geometry group. */

    struct osh_transport_scoring_page_runtime *pages; /* Pages sharing this geometry. */
    size_t npages;                                    /* Length of @ref pages. */

    size_t nbins;             /* Dense number of geometry bins. */
    unsigned int geometry_id; /* Future dense geometry enum. */
};

/**
 * @brief Transport-facing compiled scorer workspace.
 *
 * @details
 * This is the scorer analogue of `osh_transport_material_runtime`: raw parsed
 * `detect.dat` sections are compiled into dense geometry groups with attached
 * pages so the hot transport kernel can score with minimal branching.
 */
struct osh_transport_scoring_runtime {
    struct osh_transport_scoring_geometry_runtime *geometries; /* Shared scoring geometries. */
    size_t ngeometries;                                        /* Length of @ref geometries. */
};

#ifdef __cplusplus
}
#endif

#endif /* OSH_TRANSPORT_SCORING_RUNTIME_H */
