#ifndef OSH_SCORING_H
#define OSH_SCORING_H

#include <stddef.h>

#include "common/osh_rc.h"

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
 * @brief Parse a `detect.dat` file into a raw scoring workspace.
 *
 * @param[in]  path    Path to `detect.dat`.
 * @param[out] ws_out  Receives a newly allocated workspace on success.
 *
 * @returns OSH_OK on success, or a parse/allocation error code.
 */
enum osh_status osh_scoring_setup_from_path(char const *path, struct osh_scoring_workspace **ws_out);

/**
 * @brief Free a scoring workspace allocated by osh_scoring_setup_from_path().
 */
void osh_scoring_workspace_free(struct osh_scoring_workspace *ws);

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
    char *quantity;      /* Quantity keyword, e.g. "DOSE", "FLUENCE". */
    char **filter_names; /* Referenced filter names (owned). */
    size_t nfilter_names;
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
 *   Mesh  — three axes: X, Y, Z.
 *   Cyl   — two axes: R, Z.
 *   Zone  — zero axes; zone range is in @ref zone_start / @ref zone_stop.
 *
 * Rotation (theta, phi) is optional; @ref has_rotation is set when present.
 */
struct osh_scoring_geometry_def {
    char *kind;                        /* Geometry type keyword. */
    char *name;                        /* User-visible geometry name. */
    struct osh_scoring_axis_def *axes; /* Axis definitions (owned). */
    size_t naxes;                      /* Number of entries in axes[]. */
    double rot_theta_deg;              /* Rotation polar angle [deg]. */
    double rot_phi_deg;                /* Rotation azimuth angle [deg]. */
    int zone_start;                    /* First zone (Zone geometry only). */
    int zone_stop;                     /* Last zone  (Zone geometry only). */
    unsigned char has_rotation;        /* Whether rot_* fields were set. */
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
    char label[4]; /* "X", "Y", "Z", "R" — null-terminated. */
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
    double rescale;          /* Optional multiplicative output rescale. */
    double offset;           /* Optional additive output offset. */
    double site_diameter_um; /* Optional site diameter [um]. */
    double density_g_cm3;    /* Optional density override [g/cm^3]. */
    size_t npart;            /* Optional particle-count cap. */
    int medium;              /* Optional medium override. */
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

#endif /* OSH_SCORING_H */
