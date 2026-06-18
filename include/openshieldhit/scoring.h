#ifndef OPENSHIELDHIT_SCORING_H
#define OPENSHIELDHIT_SCORING_H

#include <stddef.h>
#include <stdint.h>

#include "openshieldhit/status.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---- Public Forward Declarations ----------------------------------------- */

struct osh_scoring_workspace;
struct osh_scoring_output_def;
struct osh_scoring_page_def;
struct osh_scoring_geometry_def;
struct osh_scoring_axis_def;
struct osh_scoring_settings_def;
struct osh_scoring_filter_def;
struct osh_scoring_filter_rule;

/* ---- Public API ---------------------------------------------------------- */

/**
 * @brief Allocate an empty scoring workspace.
 *
 * @param[out] ws_out  Receives the allocated workspace on success.
 *
 * @returns OSH_OK on success, OSH_ENOMEM on allocation failure.
 */
enum osh_status osh_scoring_workspace_create(struct osh_scoring_workspace **ws_out);

/**
 * @brief Free a scoring workspace and all owned parsed scoring definitions.
 *
 * @param[in] ws  Workspace to release. Safe to call with NULL.
 *
 * @returns OSH_OK.
 */
enum osh_status osh_scoring_workspace_free(struct osh_scoring_workspace *ws);

/**
 * @brief Find a parsed filter by name.  Returns NULL if not found.
 */
struct osh_scoring_filter_def const *osh_scoring_filter_by_name(struct osh_scoring_workspace const *ws,
                                                                char const *name);

/**
 * @brief Find a parsed settings section by name.  Returns NULL if not found.
 */
struct osh_scoring_settings_def const *osh_scoring_settings_by_name(struct osh_scoring_workspace const *ws,
                                                                    char const *name);

/**
 * @brief Find a parsed geometry by name.  Returns NULL if not found.
 */
struct osh_scoring_geometry_def const *osh_scoring_geometry_by_name(struct osh_scoring_workspace const *ws,
                                                                    char const *name);

/**
 * @brief Find a parsed output by file name.  Returns NULL if not found.
 */
struct osh_scoring_output_def const *osh_scoring_output_by_filename(struct osh_scoring_workspace const *ws,
                                                                    char const *filename);

/* ---- Memory estimate ----------------------------------------------------- */

/** Fixed capacity for the geometry-name hint in @ref osh_scoring_mem_estimate. */
#define OSH_SCORING_MEM_NAME_CAP 64u

/**
 * @brief Predicted memory footprint of the scoring accumulators for a run.
 *
 * @details
 * Computed from the parsed (cold) scoring configuration *before* any buffers
 * are allocated, so callers can detect an out-of-memory situation up front
 * rather than crashing mid-allocation.  It accounts only for the scoring
 * accumulator arrays (the dominant, configuration-driven allocation); the much
 * smaller geometry/material/particle-pool memory is not included.
 *
 * The figure mirrors exactly what @c osh_scoring_compile() will later allocate:
 * for every scored page, `bins × sizeof(double)` for the primary accumulator,
 * doubled for "average" quantities (LET/Qeff) that also keep a weight
 * accumulator.  See @ref osh_scoring_estimate_memory.
 */
struct osh_scoring_mem_estimate {
    uint64_t accum_bytes;                            /**< Total bytes across all scoring accumulator arrays. */
    size_t npages;                                   /**< Number of scored pages (geometry × quantity). */
    uint64_t largest_page_bytes;                     /**< Bytes of the single largest scored page. */
    char largest_geometry[OSH_SCORING_MEM_NAME_CAP]; /**< Geometry name of that page (for messages). */
};

/**
 * @brief Estimate the scoring accumulator memory for a parsed configuration.
 *
 * @details
 * Pure, allocation-free sizing pass over the parsed workspace.  Intended to be
 * called after the workspace is fully populated (including any app-side voxel
 * grid setup) but before @c osh_simulation_create(), so the caller can apply a
 * memory budget and refuse a run that would otherwise exhaust RAM.  The result
 * matches the buffers @c osh_scoring_compile() allocates, because it reuses the
 * same per-geometry bin count and the same "uses a weight accumulator" rule.
 *
 * Pages whose geometry cannot be resolved or whose bin count is zero contribute
 * nothing (they are rejected or harmless at compile time).
 *
 * @param[out] out  Receives the estimate (fields are reset by the function).
 * @returns OSH_OK on success, OSH_EINVAL if @p ws or @p out is NULL.
 */
enum osh_status osh_scoring_estimate_memory(struct osh_scoring_workspace const *ws,
                                            struct osh_scoring_mem_estimate *out);

/* ---- Top-level workspace ------------------------------------------------- */

/**
 * @brief Top-level parsed scoring configuration.
 *
 * @details
 * Raw parsed form of `detect.dat`.  Name resolution and runtime compilation
 * happen in a separate finalize step.
 */
struct osh_scoring_workspace {
    char *fname;
    struct osh_scoring_filter_def *filters;
    struct osh_scoring_settings_def *settings;
    struct osh_scoring_geometry_def *geometries;
    struct osh_scoring_output_def *outputs;
    size_t nfilters;
    size_t nsettings;
    size_t ngeometries;
    size_t noutputs;
};

/* ---- Output (estimator + pages) ----------------------------------------- */

/**
 * @brief Parsed output definition from `detect.dat`.
 *
 * @details
 * Each Output block maps to one output file.  `fileformat` is the format
 * keyword as written (e.g. "BDO", "TEXT"); NULL means unset (default BDO).
 */
struct osh_scoring_output_def {
    char *filename;      /* Output file name. */
    char *geometry_name; /* Referenced geometry name. */
    char *fileformat;    /* Optional format keyword. */
    struct osh_scoring_page_def *pages;
    size_t npages;
};

/**
 * @brief One scored quantity/page attached to an output.
 */
struct osh_scoring_page_def {
    char *quantity;      /* Quantity keyword, e.g. "DOSE", "FLUENCE" (owned). */
    char **filter_names; /* Referenced filter names (owned). */
    size_t nfilter_names;
    /* Differential axis — all zero/NULL when no differential scoring. */
    size_t diff_nbins;   /* > 0 activates differential mode. */
    double diff_lo;      /* Lower bound of the differential axis. */
    double diff_hi;      /* Upper bound of the differential axis. */
    int diff_log;        /* 0 = linear binning, 1 = logarithmic binning. */
    char *diff_kind_str; /* Axis type keyword: "ekin", "let", "qeff", etc. (owned). */
    /* Second differential axis — all zero/NULL when unused (requires diff_nbins > 0). */
    size_t diff2_nbins; /* > 0 activates double-differential mode. */
    double diff2_lo;
    double diff2_hi;
    int diff2_log;
    char *diff2_kind_str;       /* Axis type keyword (owned). */
    char *diff2_kind_sset_name; /* Optional Settings name for diff2 axis SP override (owned). */
    /* Optional Settings name that overrides the stopping-power medium used for
     * the diff1 LET/QEFF axis — written as "Diff1Type DEDX in_Si".
     * NULL means use the transport medium (default). */
    char *diff_kind_sset_name; /* Optional Settings name for diff1 axis SP override (owned). */
};

/* ---- Geometry ------------------------------------------------------------ */

/**
 * @brief Parsed scoring geometry definition from `detect.dat`.
 *
 * @details
 * The `kind` field holds the geometry type keyword as written in the file
 * (e.g. "Mesh", "Cyl", "Zone").  The axes array holds the parsed axis
 * definitions; how many axes are present depends on the geometry type:
 *
 *   Mesh  - three axes: X, Y, Z.
 *   Cyl   - two axes: R, Z.
 *   Zone  - zero axes; zone range is in @ref zone_start / @ref zone_stop.
 *
 * When @ref has_rotation is set, @ref t maps universe coordinates to the
 * geometry's local frame (same 4×4 row-major layout as the geometry body
 * transform).  Axis bounds are expressed in local coordinates.
 */
struct osh_scoring_geometry_def {
    char *kind;                        /* Geometry type keyword. */
    char *name;                        /* User-visible geometry name. */
    struct osh_scoring_axis_def *axes; /* Axis definitions (owned). */
    size_t naxes;                      /* Number of entries in axes[]. */
    double t[16];                      /* Universe→local affine transform (row-major 4×4). */
    int zone_start;                    /* First zone (Zone geometry only). */
    int zone_stop;                     /* Last zone  (Zone geometry only). */
    unsigned char has_rotation;        /* When set, t[] is valid and axes are in local frame. */
    /* ---- Voxel geometry fields (DicomCT / DicomRTDOSE) ------------------- */
    char *vox_rtdose_path; /* DicomRTDOSE only: path to RTDOSE DICOM file; read by app, never library. */
    char *vox_body_name;   /* Future multi-CT: explicit CT body reference. NULL = auto-detect single body. */
    double vox_origin[3];  /* Voxel-grid corner [cm] in universe frame; set by app. */
    double vox_spacing[3]; /* Voxel spacing [cm]; set by app. */
    size_t vox_nx;         /* Grid dimension X; set by app. */
    size_t vox_ny;         /* Grid dimension Y; set by app. */
    size_t vox_nz;         /* Grid dimension Z; set by app. */
};

/**
 * @brief One parsed axis of a scoring geometry.
 *
 * @details
 * For Mesh geometries three axes (X, Y, Z) are expected.
 * For Cyl geometries two (R, Z).
 * For Zone geometries none (the zone list is in osh_scoring_geometry_def.zones).
 *
 * A negative nbins means "inherit from output" (old -1 convention), kept for
 * compatibility; the runtime finalize step resolves this.
 */
struct osh_scoring_axis_def {
    char label[4]; /* "X", "Y", "Z", "R" - null-terminated. */
    double lo;     /* Lower bound [cm]. */
    double hi;     /* Upper bound [cm]. */
    int nbins;     /* Number of bins; negative means unset. */
};

/* ---- Settings ------------------------------------------------------------ */

/**
 * @brief Parsed settings definition from `detect.dat`.
 */
struct osh_scoring_settings_def {
    char *name;              /* User-visible settings name. */
    char *material_name;     /* Material name from "Material <name>"; NULL if not set. */
    double rescale;          /* Optional multiplicative output rescale. */
    double offset;           /* Optional additive output offset. */
    double site_diameter_um; /* Optional site diameter [um]. */
    double density_g_cm3;    /* Optional density override [g/cm^3]. */
    size_t npart;            /* Optional particle-count cap. */
    int medium;              /* Optional medium override (dense material index). */
    int nkmedium;            /* Optional neutron-kerma medium override. */
    char has_rescale;
    char has_offset;
    char has_site_diameter_um;
    char has_density_g_cm3;
    char has_npart;
    char has_medium;
    char has_nkmedium;
};

/* ---- Filter -------------------------------------------------------------- */

/**
 * @brief Parsed filter definition from `detect.dat`.
 */
struct osh_scoring_filter_def {
    char *name;                            /* User-visible filter name. */
    struct osh_scoring_filter_rule *rules; /* Array of filter rules. */
    size_t nrules;                         /* Number of rules. */
};

/**
 * @brief One rule within a filter definition.
 *
 * @details
 * Raw text form as parsed: field name, operator string, value string.
 * The finalize step resolves field and operator to integer codes.
 *
 * Supported field names (case-insensitive): Z, A, E, GEN, ID, ...
 * Supported operators: =, <, >, <=, >=, !=
 */
struct osh_scoring_filter_rule {
    char field[16]; /* Field keyword, e.g. "Z", "A", "E", "GEN". */
    char op[4];     /* Operator string, e.g. "=", ">", "<=". */
    double value;   /* Numeric right-hand side. */
};

#ifdef __cplusplus
}
#endif

#endif /* OPENSHIELDHIT_SCORING_H */
