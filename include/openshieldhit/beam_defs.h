#ifndef OPENSHIELDHIT_BEAM_DEFS_H
#define OPENSHIELDHIT_BEAM_DEFS_H

#define OSH_BEAM_STRAGG_OFF 0
#define OSH_BEAM_STRAGG_GAUSS 1
#define OSH_BEAM_STRAGG_VAVILOV 2
#define OSH_BEAM_STRAGG_URBAN 3

#define OSH_BEAM_MSCAT_OFF 0
#define OSH_BEAM_MSCAT_GAUSS 1
#define OSH_BEAM_MSCAT_MOLIERE 2
#define OSH_BEAM_MSCAT_WENTZEL 3

#define OSH_BEAM_MODE_SPOTS 0
#define OSH_BEAM_MODE_SOBP 1
#define OSH_BEAM_MODE_PHSP 2

#define OSH_BEAM_SHAPE_PENCIL 0
#define OSH_BEAM_SHAPE_GAUSSIAN 1
#define OSH_BEAM_SHAPE_SQUARE 2
#define OSH_BEAM_SHAPE_CIRCULAR 3
#define OSH_BEAM_SHAPE_INVALID 255

/* Minimum primary-beam kinetic energy [MeV/nucleon].  Deliberately higher than
 * the transport tracking cutoff (OSH_MATERIAL_RUNTIME_EMIN = 0.025 MeV/u, the
 * LOADDEDX table floor at which in-flight particles are killed): that cutoff
 * governs when a tracked particle stops, whereas this is the lowest energy a
 * *primary* may be launched at.  Beams below 0.1 MeV/u have sub-millimetre range
 * and no practical use, so they are rejected at parse as likely user error. */
#define OSH_BEAM_TMIN 0.1

static char const *const osh_beam_stragg_names[] = {"OFF", "GAUSSIAN", "VAVILOV", "URBAN"};
static char const *const osh_beam_mscat_names[] = {"OFF", "GAUSSIAN", "MOLIERE", "WENTZEL"};
static char const *const osh_beam_mode_names[] = {"SPOTS", "SOBP", "PHASESPACE"};
static char const *const osh_beam_shape_names[] = {"PENCIL", "GAUSSIAN", "SQUARE", "CIRCULAR"};

#endif /* OPENSHIELDHIT_BEAM_DEFS_H */
