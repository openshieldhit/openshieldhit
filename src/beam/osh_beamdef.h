#ifndef OSH_BEAMDEF_H
#define OSH_BEAMDEF_H

/* ---- Beam physics-model constants ----------------------------------------
 *
 * Kept in a separate header so that output/metadata code can include just
 * this file without pulling in the full beam workspace structs. */

#define OSH_BEAM_STRAGG_OFF 0     /* no energy straggling */
#define OSH_BEAM_STRAGG_GAUSS 1   /* Gaussian straggling */
#define OSH_BEAM_STRAGG_VAVILOV 2 /* Vavilov straggling */

#define OSH_BEAM_MSCAT_OFF 0     /* no multiple scattering */
#define OSH_BEAM_MSCAT_GAUSS 1   /* Gaussian (Highland) model */
#define OSH_BEAM_MSCAT_MOLIERE 2 /* Moliere model */

#define OSH_BEAM_MODE_SPOTS 0 /* pencil/scanning spot list */
#define OSH_BEAM_MODE_SOBP 1  /* spread-out Bragg peak from external file */
#define OSH_BEAM_MODE_PHSP 2  /* phase-space source (MCPL) */

#define OSH_BEAM_SHAPE_PENCIL 0   /* zero-width pencil beam */
#define OSH_BEAM_SHAPE_GAUSSIAN 1 /* Gaussian lateral profile */
#define OSH_BEAM_SHAPE_SQUARE 2   /* uniform square profile */
#define OSH_BEAM_SHAPE_CIRCULAR 3 /* uniform circular/annular profile */
#define OSH_BEAM_SHAPE_INVALID 255

/* Minimum beam energy at which transport is meaningful [MeV total]. */
#define OSH_BEAM_TMIN 0.1

/* ---- Human-readable name arrays ------------------------------------------
 *
 * Indexed directly by the OSH_BEAM_* constants above.
 * Declared static so each translation unit gets its own copy — acceptable
 * for these small arrays.  Use these wherever beam mode/shape names appear
 * in output or metadata files so the strings stay in sync with the defines. */

static char const *const osh_beam_stragg_names[] = {"OFF", "GAUSSIAN", "VAVILOV"};

static char const *const osh_beam_mscat_names[] = {"OFF", "GAUSSIAN", "MOLIERE"};

static char const *const osh_beam_mode_names[] = {"SPOTS", "SOBP", "PHASESPACE"};

static char const *const osh_beam_shape_names[] = {"PENCIL", "GAUSSIAN", "SQUARE", "CIRCULAR"};

#endif /* OSH_BEAMDEF_H */
