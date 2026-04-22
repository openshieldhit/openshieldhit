#ifndef OSH_APP_OSH_GEOMETRY_PARSE_KEYS_H
#define OSH_APP_OSH_GEOMETRY_PARSE_KEYS_H

/**
 * @file osh_geometry_parse_keys.h
 * @brief OpenShieldHIT `geo.dat` file-syntax keys used by the app parser.
 *
 * @details
 * These are app-private parsing tokens. They intentionally do not live in
 * `osh_core`, because file syntax belongs to `src/apps/osh/`, not to the
 * public geometry API or the GEMCA prepare/runtime layers.
 */

#define OSH_GEO_KEY_ASSIGNMA "assignma"   /* Legacy ASSIGNMAT alias. */
#define OSH_GEO_KEY_ASSIGNMAT "assignmat" /* Material assignment card. */
#define OSH_GEO_KEY_END "end"             /* Section delimiter card. */

/* List of known body-card keys (geo.dat body section). */
#define OSH_GEO_KEY_ARB "arb" /* Arbitrary convex polyhedron body card. */
#define OSH_GEO_KEY_BOX "box" /* Box body card. */
#define OSH_GEO_KEY_CPY "cpy" /* Body copy helper card. */
#define OSH_GEO_KEY_DCM "dcm" /* CT DICOM-backed voxel convenience card. */
#define OSH_GEO_KEY_ELL "ell" /* Ellipsoid body card. */

#define OSH_GEO_KEY_MOV "mov" /* Body move helper card. */
#define OSH_GEO_KEY_PLA "pla" /* Generic plane card. */
#define OSH_GEO_KEY_RCC "rcc" /* Right circular cylinder body card. */
#define OSH_GEO_KEY_REC "rec" /* Right elliptical cylinder body card. */
#define OSH_GEO_KEY_ROT "rot" /* Body rotate helper card. */
#define OSH_GEO_KEY_RPP "rpp" /* Axis-aligned RPP body card. */
#define OSH_GEO_KEY_SPH "sph" /* Sphere body card. */
#define OSH_GEO_KEY_TRC "trc" /* Truncated cone body card. */
#define OSH_GEO_KEY_VOX "vox" /* Voxel body card. */
#define OSH_GEO_KEY_WED "wed" /* Wedge body card. */
#define OSH_GEO_KEY_XYP "xyp" /* Plane normal to Z (XY plane). */
#define OSH_GEO_KEY_XZP "xzp" /* Plane normal to Y (XZ plane). */
#define OSH_GEO_KEY_YZP "yzp" /* Plane normal to X (YZ plane). */

#endif /* OSH_APP_OSH_GEOMETRY_PARSE_KEYS_H */
