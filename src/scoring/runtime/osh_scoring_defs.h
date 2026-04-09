#ifndef OSH_SCORING_DEFS_H
#define OSH_SCORING_DEFS_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Shared scoring enumeration types.
 *
 * This header is included by both geometry and output/page runtime headers
 * to avoid any circular dependency between those two layers.  It must not
 * include any other scoring headers.
 *
 * For the legacy equivalent see _temp_shieldhit/scoring/sh_scoredef.h,
 * which uses #define constants instead of enums.
 */

/* ---- Geometry type ------------------------------------------------------ */

/**
 * @brief Compiled geometry kind for a scoring region.
 *
 * @details
 * Resolved from the keyword string in detect.dat during osh_scoring_prepare().
 * The integer codes are used in the hot-path geometry dispatch so string
 * comparisons are not needed at score time.
 *
 * Mesh and Cyl geometries use axis-aligned bin lookup (see
 * osh_scoring_axis_runtime).  Zone geometry delegates zone membership to
 * the GEMCA engine.  Voxel geometry uses the raytrace module
 * (common/raytrace/) for step-length decomposition.
 */
enum osh_scoring_geo_kind {
    OSH_SCORING_GEO_UNKNOWN = 0, /* not yet resolved — should not reach score_step */
    OSH_SCORING_GEO_MESH    = 1, /* Cartesian (X,Y,Z) mesh */
    OSH_SCORING_GEO_CYL     = 2, /* cylindrical (R,Z) mesh */
    OSH_SCORING_GEO_ZONE    = 3, /* GEMCA zone scoring */
    OSH_SCORING_GEO_VOXEL   = 4, /* CT voxel grid (uses raytrace) */
    OSH_SCORING_GEO_ALL     = 5  /* whole simulation universe */
};

/* ---- Postprocessing mode ------------------------------------------------ */

/**
 * @brief How accumulated page data is combined across simulation runs.
 *
 * @details
 * After N simulation runs (instances), each contributing x_j scored with
 * I_j primary particles, the final result X is derived as follows:
 *
 *   NONE   — X = x_0
 *              Used for geometry-map scorers (GEOMAP) where only one pass
 *              is meaningful.
 *
 *   SUM    — X = sum_j x_j
 *              Raw sum across runs.  Used for COUNT, MCPL particle count.
 *
 *   NORM   — X = (sum_j x_j) / (sum_j I_j)
 *              Normalised per primary particle.  Used for DOSE, FLUENCE,
 *              NORMCOUNT, and most physical quantity scorers.
 *
 *   AVER   — X = (sum_j x_j * I_j) / (sum_j I_j)
 *              Weighted average across runs.  Used for LET, averaged
 *              quantities where each run already produces a per-event mean.
 *
 *   APPEND — X = [x_0, x_1, ..., x_j]
 *              Sequential concatenation.  Used for MCPL phase-space output
 *              where each particle record is appended rather than binned.
 */
enum osh_scoring_postproc {
    OSH_SCORING_POSTPROC_NONE   = 0,
    OSH_SCORING_POSTPROC_SUM    = 1,
    OSH_SCORING_POSTPROC_NORM   = 2,
    OSH_SCORING_POSTPROC_AVER   = 3,
    OSH_SCORING_POSTPROC_APPEND = 4
};

/* ---- Scored quantity (detector type) ------------------------------------ */

/**
 * @brief The physical quantity scored into a page's data array.
 *
 * @details
 * Each page accumulates exactly one scored quantity.  The kind determines
 * which algorithm runs in score_step/score_point and what physical units
 * the data array holds before postprocessing.
 *
 * This list will grow as detectors are implemented.  For the full legacy
 * list (63 entries) see _temp_shieldhit/scoring/sh_scoredef.h SH_SDET_*.
 *
 * Naming conventions
 * ------------------
 * D-prefix  : dose-averaged quantity  (e.g. DLET = dose-averaged LET)
 * T-prefix  : track-averaged quantity (e.g. TLET = track-averaged LET)
 * No prefix : fluence-like or count-like (un-averaged)
 *
 * Units after postprocessing (NORM mode unless noted)
 * ---------------------------------------------------
 * DOSE    [MeV/g]       multiply by 1.602e-10 to get [Gy]
 * FLUENCE [1/cm^2]
 * LET     [MeV/cm]      or [keV/um] depending on settings rescale
 * COUNT   [1]           raw event count (no normalisation)
 * NKERMA  [MeV/g]       neutron kerma, same unit as DOSE
 */
enum osh_scoring_score_kind {

    OSH_SCORING_SCORE_UNKNOWN  =  0, /* not yet resolved */

    /* --- Basic transport quantities ------------------------------------- */

    OSH_SCORING_SCORE_ENERGY   =  1, /* mean energy deposited in voxel [MeV] */
    OSH_SCORING_SCORE_FLUENCE  =  2, /* particle fluence (Chilton / ICRU definition) [1/cm^2] */
    OSH_SCORING_SCORE_DOSE     =  3, /* absorbed dose [MeV/g]; use rescale for [Gy] */
    OSH_SCORING_SCORE_LETFLU   =  4, /* LET * fluence [MeV/cm * 1/cm^2] — numerator for TLET */

    /* --- LET-averaged quantities ---------------------------------------- */

    OSH_SCORING_SCORE_DLET     =  5, /* dose-averaged LET  (Cortez algorithm C) [MeV/cm] */
    OSH_SCORING_SCORE_TLET     =  6, /* track-averaged LET (Cortez algorithm C) [MeV/cm] */

    /* --- Counters ------------------------------------------------------- */

    OSH_SCORING_SCORE_NORMCOUNT =  7, /* event counter, normalised per primary [1] */
    OSH_SCORING_SCORE_COUNT     =  8, /* event counter, not normalised [1] */

    /* --- Neutron quantities --------------------------------------------- */

    OSH_SCORING_SCORE_NKERMA   =  9, /* neutron kerma [MeV/g] */

    /* --- Detector models ------------------------------------------------ */

    OSH_SCORING_SCORE_ALANINE  = 10, /* alanine detector response (Bassler 2008) [1] */

    /* --- Phase-space output --------------------------------------------- */

    OSH_SCORING_SCORE_MCPL     = 11  /* MC particle list (phase-space append mode) */
};

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_DEFS_H */
