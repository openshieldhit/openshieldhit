#ifndef OSH_TRANSPORT_SCORING_RUNTIME_H
#define OSH_TRANSPORT_SCORING_RUNTIME_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Compiled filter and settings metadata -------------------------------- */

struct osh_transport_scoring_filter_rule {
    char field[16];
    char op[4];
    double value;
};

struct osh_transport_scoring_filter_runtime {
    char *name;
    struct osh_transport_scoring_filter_rule *rules;
    size_t nrules;
};

struct osh_transport_scoring_settings_runtime {
    char *name;
    double rescale;
    double offset;
    double site_diameter_um;
    double density_g_cm3;
    size_t npart;
    int medium;
    int nkmedium;
    unsigned char has_rescale;
    unsigned char has_offset;
    unsigned char has_site_diameter_um;
    unsigned char has_density_g_cm3;
    unsigned char has_npart;
    unsigned char has_medium;
    unsigned char has_nkmedium;
};

/* ---- Compiled geometry ---------------------------------------------------- */

struct osh_transport_scoring_axis_runtime {
    char label[4];
    double lo;
    double hi;
    int nbins;
};

struct osh_transport_scoring_geometry_runtime {
    char *kind;
    char *name;
    struct osh_transport_scoring_axis_runtime *axes;
    size_t naxes;
    double rot_theta_deg;
    double rot_phi_deg;
    int zone_start;
    int zone_stop;
    size_t nbins;
    size_t first_page;
    size_t npages;
    unsigned char has_rotation;
};

/* ---- Page/runtime scoring buffers ---------------------------------------- */

struct osh_transport_scoring_page_filter_ref {
    size_t filter_idx;
};

struct osh_transport_scoring_page_settings_ref {
    size_t settings_idx;
};

struct osh_transport_scoring_page_runtime {
    char *quantity;
    struct osh_transport_scoring_page_filter_ref *filters;
    struct osh_transport_scoring_page_settings_ref *settings;
    double *data;
    double *data_var;
    double *data2;
    double *data2_var;
    size_t output_idx;
    size_t geometry_idx;
    size_t nfilters;
    size_t nsettings;
    size_t len;
    unsigned char has_data2;
    unsigned char variance;
    unsigned char divide;
    unsigned char postproc;
};

/* ---- Output/runtime file grouping ---------------------------------------- */

struct osh_transport_scoring_output_runtime {
    char *filename;
    char *fileformat;
    size_t geometry_idx;
    size_t *page_indices;
    size_t npages;
};

/* ---- Top-level compiled runtime ------------------------------------------ */

/**
 * @brief Transport-facing compiled scorer workspace.
 *
 * @details
 * The prepare step resolves raw `detect.dat` names to dense indices and
 * allocates page-local accumulators. Runtime traversal is geometry-centric:
 * each geometry owns one contiguous page span so Siddon/Jacobs traversal can
 * be executed once per geometry, then all attached pages updated in one loop.
 *
 * Output files remain a separate cold-path grouping. Each output stores the
 * page indices that should later be written into one BDO file/page sequence.
 */
struct osh_transport_scoring_runtime {
    struct osh_transport_scoring_filter_runtime *filters;
    struct osh_transport_scoring_settings_runtime *settings;
    struct osh_transport_scoring_geometry_runtime *geometries;
    struct osh_transport_scoring_page_runtime *pages;
    struct osh_transport_scoring_output_runtime *outputs;
    size_t nfilters;
    size_t nsettings;
    size_t ngeometries;
    size_t npages;
    size_t noutputs;
};

#ifdef __cplusplus
}
#endif

#endif /* OSH_TRANSPORT_SCORING_RUNTIME_H */
