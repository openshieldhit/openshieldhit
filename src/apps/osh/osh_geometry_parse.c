#include "apps/osh/osh_geometry_parse.h"

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_geometry_parse_keys.h"
#include "common/osh_diag.h"
#include "common/osh_file.h"
#include "common/osh_patient_position.h"
#include "common/osh_readline.h"
#include "common/osh_vect.h"
#include "common/osh_voxel_order.h"
#include "gemca/osh_gemca2_defines.h"
#include "openshieldhit/const.h"
#include "openshieldhit/dicom.h"
#include "openshieldhit/geometry.h"
#include "openshieldhit/geometry_defs.h"
#include "openshieldhit/status.h"

/* ---- Internal helpers ---------------------------------------------------- */

static int _rewind_oshfile(struct oshfile *shf, struct osh_diag_sink const *diag);
static enum osh_status _finalize_body(
    struct osh_geometry_workspace *ws, size_t ibody, int btype, char const *name, double const *par, int npar);
static int _parse_double_token(char const *token, double *out);
static enum osh_status _parse_dcm_body(char const *args,
                                       char const *geo_filename,
                                       int lineno,
                                       struct osh_diag_sink const *diag,
                                       char name_out[OSH_GEMCA_BODY_NAME_MAXLEN],
                                       double par_out[OSH_GEMCA_NARGS_MAX],
                                       int *npar_out,
                                       int16_t **hu_out,
                                       size_t *n_hu_out);
static int _ct_orientation_is_axial(struct osh_dicom_ct const *ct);
static int _mul_size_overflow(size_t a, size_t b, size_t *out);

/**
 * @brief Map a body-type key string to an OSH_GEOMETRY_BODY_* code.
 *
 * @details Values are identical to the internal OSH_GEMCA_BODY_* constants so
 * that osh_geometry_workspace_prepare() can use them without conversion.
 *
 * @param[in] key  Lowercase or mixed-case body type token (e.g. "sph", "rpp").
 *
 * @returns The matching OSH_GEOMETRY_BODY_* constant, or OSH_GEOMETRY_BODY_NONE
 *          if @p key is unrecognized.
 */
static int _body_type_from_key(char const *key) {
    if (strcasecmp(key, OSH_GEO_KEY_SPH) == 0) {
        return OSH_GEOMETRY_BODY_SPH;
    }
    if (strcasecmp(key, OSH_GEO_KEY_WED) == 0) {
        return OSH_GEOMETRY_BODY_WED;
    }
    if (strcasecmp(key, OSH_GEO_KEY_ARB) == 0) {
        return OSH_GEOMETRY_BODY_ARB;
    }
    if (strcasecmp(key, OSH_GEO_KEY_BOX) == 0) {
        return OSH_GEOMETRY_BODY_BOX;
    }
    if (strcasecmp(key, OSH_GEO_KEY_VOX) == 0) {
        return OSH_GEOMETRY_BODY_VOX;
    }
    if (strcasecmp(key, OSH_GEO_KEY_DCM) == 0) {
        return OSH_GEOMETRY_BODY_VOX;
    }
    if (strcasecmp(key, OSH_GEO_KEY_RPP) == 0) {
        return OSH_GEOMETRY_BODY_RPP;
    }
    if (strcasecmp(key, OSH_GEO_KEY_RCC) == 0) {
        return OSH_GEOMETRY_BODY_RCC;
    }
    if (strcasecmp(key, OSH_GEO_KEY_REC) == 0) {
        return OSH_GEOMETRY_BODY_REC;
    }
    if (strcasecmp(key, OSH_GEO_KEY_TRC) == 0) {
        return OSH_GEOMETRY_BODY_TRC;
    }
    if (strcasecmp(key, OSH_GEO_KEY_ELL) == 0) {
        return OSH_GEOMETRY_BODY_ELL;
    }
    if (strcasecmp(key, OSH_GEO_KEY_YZP) == 0) {
        return OSH_GEOMETRY_BODY_YZP;
    }
    if (strcasecmp(key, OSH_GEO_KEY_XZP) == 0) {
        return OSH_GEOMETRY_BODY_XZP;
    }
    if (strcasecmp(key, OSH_GEO_KEY_XYP) == 0) {
        return OSH_GEOMETRY_BODY_XYP;
    }
    if (strcasecmp(key, OSH_GEO_KEY_PLA) == 0) {
        return OSH_GEOMETRY_BODY_PLA;
    }
    if (strcasecmp(key, OSH_GEO_KEY_ROT) == 0) {
        return OSH_GEOMETRY_BODY_ROT;
    }
    if (strcasecmp(key, OSH_GEO_KEY_CPY) == 0) {
        return OSH_GEOMETRY_BODY_CPY;
    }
    if (strcasecmp(key, OSH_GEO_KEY_MOV) == 0) {
        return OSH_GEOMETRY_BODY_MOV;
    }
    return OSH_GEOMETRY_BODY_NONE;
}

/**
 * @brief Check whether @p key is a zone boolean-expression continuation token.
 *
 * @param[in] key  Token to test.
 *
 * @returns 1 if @p key starts with +, -, |, ( or ), 0 otherwise.
 */
static int _is_zone_continuation(char const *key) {
    return key && (key[0] == '+' || key[0] == '-' || key[0] == '|' || key[0] == '(' || key[0] == ')');
}

/**
 * @brief Count body cards in the geometry file without retaining state.
 *
 * @param[in,out] shf  Open geometry file; rewound internally.
 *
 * @returns Number of body cards found, or 0 on I/O error.
 */
static size_t _count_bodies(struct oshfile *shf, struct osh_diag_sink const *diag) {
    int lineno;
    char *line = NULL;
    char *key = NULL;
    char *args = NULL;
    size_t nbody = 0u;

    if (!_rewind_oshfile(shf, diag)) {
        return 0u;
    }

    while (osh_readline_key(shf, &line, &key, &args, &lineno) > 0) {
        if (strcasecmp(key, OSH_GEO_KEY_END) != 0 && _body_type_from_key(key) != OSH_GEOMETRY_BODY_NONE) {
            ++nbody;
        }
        free(line);
        line = NULL;
    }

    free(line);
    return nbody;
}

/**
 * @brief Count zones in the geometry file without building any cold objects.
 *
 * @param[in,out] shf  Open geometry file; rewound internally.
 *
 * @returns Number of zone header cards found, or 0 on I/O error.
 */
static size_t _count_zones(struct oshfile *shf, struct osh_diag_sink const *diag) {
    int lineno;
    char *line = NULL;
    char *key = NULL;
    char *args = NULL;
    size_t nzone = 0u;
    int in_block = 0;

    if (!_rewind_oshfile(shf, diag)) {
        return 0u;
    }

    while (osh_readline_key(shf, &line, &key, &args, &lineno) > 0) {
        if (strcasecmp(key, OSH_GEO_KEY_END) != 0) {
            if (in_block == 2 && !_is_zone_continuation(key)) {
                ++nzone;
            }
        } else {
            in_block++;
        }

        if (in_block == 0) {
            in_block = 1;
        }

        free(line);
        line = NULL;
    }

    free(line);
    OSH_DIAG_INFOF(diag, "Found %zu zones in geo.dat file", nzone);
    return nzone;
}

/**
 * @brief Append string @p b to the heap-allocated string at @p *a.
 *
 * @param[in,out] a  Pointer to the destination heap string; grown as needed.
 * @param[in]     b  NUL-terminated string to append.
 *
 * @returns OSH_OK on success, OSH_ENOMEM on allocation failure.
 */
static enum osh_status _str_append(char **a, char const *b) {
    size_t la = strlen(*a);
    size_t lb = strlen(b);
    char *p = (char *) realloc(*a, la + lb + 1u);
    if (!p) {
        return OSH_ENOMEM;
    }
    memcpy(p + la, b, lb + 1u);
    *a = p;
    return OSH_OK;
}

/**
 * @brief Find a zone by exact name in the flat cold-zone array.
 *
 * @param[in]  name     Zone name to look up.
 * @param[in]  zones    Array of cold zone definitions.
 * @param[in]  nzones   Length of @p zones.
 * @param[out] idx_out  Receives the index of the matching zone.
 *
 * @returns 1 if the name is found (fills @p *idx_out), 0 otherwise.
 */
static int
_zone_index_from_name(char const *name, struct osh_geometry_zone const *zones, size_t nzones, size_t *idx_out) {
    size_t i;
    for (i = 0u; i < nzones; ++i) {
        if (zones[i].name && strcmp(zones[i].name, name) == 0) {
            *idx_out = i;
            return 1;
        }
    }
    return 0;
}

/**
 * @brief Return the next whitespace-delimited token from @p *cursor.
 *
 * @details Advances @p *cursor past the returned token. Modifies the buffer
 * in place by inserting a NUL terminator after each token.
 *
 * @param[in,out] cursor  Pointer into the buffer; updated to point past the token.
 *
 * @returns Pointer to the start of the token, or NULL when no more tokens remain.
 */
static char *_next_token(char **cursor) {
    char *start;
    char *end;

    if (!cursor || !*cursor) {
        return NULL;
    }
    start = *cursor;
    while (*start != '\0' && isspace((unsigned char) *start)) {
        start++;
    }
    if (*start == '\0') {
        *cursor = start;
        return NULL;
    }
    end = start;
    while (*end != '\0' && !isspace((unsigned char) *end)) {
        end++;
    }
    if (*end != '\0') {
        *end = '\0';
        end++;
    }
    *cursor = end;
    return start;
}

/**
 * @brief Copy one parsed body record into the workspace body array.
 */
static enum osh_status _finalize_body(
    struct osh_geometry_workspace *ws, size_t ibody, int btype, char const *name, double const *par, int npar) {
    char *name_copy = NULL;
    double *a_copy = NULL;

    if (!ws || ibody >= ws->nbodies || !name || npar < 0 || (npar > 0 && !par)) {
        return OSH_EINVAL;
    }
    /* Allocate into temporaries first; only write to ws on full success to
     * avoid partial initialisation that would complicate caller cleanup. */
    name_copy = strdup(name);
    if (!name_copy) {
        return OSH_ENOMEM;
    }
    if (npar > 0) {
        a_copy = (double *) calloc((size_t) npar, sizeof(double));
        if (!a_copy) {
            free(name_copy);
            return OSH_ENOMEM;
        }
        memcpy(a_copy, par, (size_t) npar * sizeof(double));
    }
    ws->bodies[ibody].type = btype;
    ws->bodies[ibody].name = name_copy;
    ws->bodies[ibody].a = a_copy;
    ws->bodies[ibody].na = npar;
    return OSH_OK;
}

/**
 * @brief Parse one token as a finite double value.
 */
static int _parse_double_token(char const *token, double *out) {
    char *end = NULL;
    double v;
    if (!token || !out) {
        return 0;
    }
    errno = 0;
    v = strtod(token, &end);
    if (errno != 0 || end == token || *end != '\0' || !isfinite(v)) {
        return 0;
    }
    *out = v;
    return 1;
}

/**
 * @brief Check whether CT ImageOrientationPatient describes an axial stack.
 *
 * @details
 * Current DCM->VOX placement supports axial CT only. Axial means row/column
 * direction cosines lie in the XY plane and their cross-product normal points
 * along +/-Z. This still allows in-plane rotation around Z.
 *
 * Future work: remove this guard once full orientation-cosine placement is
 * implemented for arbitrary oblique datasets.
 */
static int _ct_orientation_is_axial(struct osh_dicom_ct const *ct) {
    double const eps = 1.0e-3;
    double const *r;
    double const *c;
    double row_norm_sq;
    double col_norm_sq;
    double row_col_dot;
    double n[3];

    if (!ct) {
        return 0;
    }

    r = ct->row_cosine;
    c = ct->col_cosine;
    row_norm_sq = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
    col_norm_sq = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
    row_col_dot = r[0] * c[0] + r[1] * c[1] + r[2] * c[2];
    n[0] = r[1] * c[2] - r[2] * c[1];
    n[1] = r[2] * c[0] - r[0] * c[2];
    n[2] = r[0] * c[1] - r[1] * c[0];

    if (fabs(row_norm_sq - 1.0) > eps || fabs(col_norm_sq - 1.0) > eps || fabs(row_col_dot) > eps) {
        return 0;
    }
    if (fabs(r[2]) > eps || fabs(c[2]) > eps) {
        return 0;
    }
    if (fabs(n[0]) > eps || fabs(n[1]) > eps || fabs(fabs(n[2]) - 1.0) > eps) {
        return 0;
    }

    return 1;
}

static int _mul_size_overflow(size_t a, size_t b, size_t *out) {
    if (!out) {
        return 1;
    }
    if (a != 0u && b > SIZE_MAX / a) {
        return 1;
    }
    *out = a * b;
    return 0;
}

/**
 * @brief Parse one DCM card and convert CT metadata into VOX raw parameters.
 *
 * @details
 * Expected card payload (after key):
 * `name ct_dir patient_pos gantry_deg couch_deg iso_x_mm iso_y_mm iso_z_mm`
 *
 * Where @c patient_pos is an IEC 61217 patient position string such as
 * "HFS", "HFP", "FFS", "FFP", "HFDL", "HFDR", "FFDL", or "FFDR".
 * iso_{x,y,z}_mm is the treatment isocenter in DICOM LPS patient coordinates [mm].
 *
 * Generated VOX parameters (18 total):
 *   0..2   x0,y0,z0  voxel-grid corner [cm] in local voxel coordinates (always zero)
 *   3..5   dx,dy,dz  voxel spacing [cm]
 *   6..8   nx,ny,nz  voxel counts
 *   9..10  gantry,couch [deg]
 *   11..13 t_corner[0..2] [cm]: universe position of the CT voxel-grid corner,
 *          computed as -(tb^T * iso_ct_local), where iso_ct_local is the isocenter
 *          expressed in CT-local frame (cm from corner).
 *   14..16 -ct.origin[0..2]/10 [cm]: negative DICOM origin offset, used by the
 *          RTDOSE/DicomCT scoring geometry to convert absolute DICOM coordinates
 *          to CT-local coordinates.  Independent of isocenter and rotation.
 *   17     patient position code (enum osh_patient_position cast to double)
 *
 * On success this function also transfers ownership of the DICOM HU pixel
 * buffer to @p hu_out (workspace-level ownership).
 *
 * Coordinate-chain (IEC 61217 universe):
 * - The simulation universe follows IEC 61217: X=patient-left, Y=cranial
 *   (head-first), Z=anterior/nozzle direction at gantry 0.
 * - The CT reader gives geometry in DICOM LPS patient coordinates.
 * - The patient-position base rotation tb maps universe coords to DICOM coords:
 *   DICOM[i] = tb[i] . p_universe.  Couch and gantry rotations are then applied
 *   to each row of tb (couch around universe Z, gantry around universe Y).
 * - iso_ct_local[j] = (iso_mm[j] - ct.origin[j]) / 10 + spacing[j] / 2
 * - t_corner[k] = -sum_j(tb[j][k] * iso_ct_local[j])  (maps iso to universe origin)
 * - Current limitation: only axial CT orientation is accepted. Non-axial
 *   datasets are rejected with a parse error.
 *
 * See also docs/voxel_coordinates.md (single-source convention note).
 */
static enum osh_status _parse_dcm_body(char const *args,
                                       char const *geo_filename,
                                       int lineno,
                                       struct osh_diag_sink const *diag,
                                       char name_out[OSH_GEMCA_BODY_NAME_MAXLEN],
                                       double par_out[OSH_GEMCA_NARGS_MAX],
                                       int *npar_out,
                                       int16_t **hu_out,
                                       size_t *n_hu_out) {
    char *work = NULL;
    char *cursor = NULL;
    char *tok = NULL;
    char *ct_dir_rel = NULL;
    char *geo_dir = NULL;
    char *ct_dir_abs = NULL;
    struct osh_dicom_ct ct;
    enum osh_status rc = OSH_EPARSE;
    enum osh_patient_position patient_pos;
    double tb[3][3]; /* base rotation: maps universe -> DICOM LPS (row i = DICOM axis i in universe) */
    double gantry_deg;
    double couch_deg;
    double gantry_rad;
    double couch_rad;
    double iso_mm[3];       /* treatment isocenter in DICOM patient coordinate system (LPS, mm) */
    double iso_ct_local[3]; /* isocenter in CT-local frame (cm, relative to CT corner) */
    int j;
    int k;

    if (!args || !name_out || !par_out || !npar_out) {
        return OSH_EINVAL;
    }
    if (!hu_out || !n_hu_out) {
        return OSH_EINVAL;
    }
    *hu_out = NULL;
    *n_hu_out = 0u;

    memset(&ct, 0, sizeof(ct));
    work = strdup(args);
    if (!work) {
        return OSH_ENOMEM;
    }
    cursor = work;

    tok = _next_token(&cursor);
    if (!tok) {
        OSH_DIAG_ERRORF(diag, "%s line %d: DCM card requires a body name", geo_filename, lineno);
        goto done;
    }
    if (strlen(tok) >= OSH_GEMCA_BODY_NAME_MAXLEN) {
        OSH_DIAG_ERRORF(diag,
                        "%s line %d: DCM body name too long (max %d chars)",
                        geo_filename,
                        lineno,
                        OSH_GEMCA_BODY_NAME_MAXLEN - 1);
        goto done;
    }
    memcpy(name_out, tok, strlen(tok) + 1u);

    ct_dir_rel = _next_token(&cursor);
    if (!ct_dir_rel) {
        OSH_DIAG_ERRORF(diag, "%s line %d: DCM card requires a CT directory path", geo_filename, lineno);
        goto done;
    }

    /* Patient position token (IEC 61217): "HFS", "HFP", "FFS", "FFP",
     * "HFDL", "HFDR", "FFDL", or "FFDR".
     * This determines the rotation from DICOM LPS patient coords to the IEC
     * universe frame.  Without it, the voxel grid would be placed with the
     * wrong axis orientation (e.g. beam entering from the cranial end instead
     * of the anterior face for a standard brain treatment). */
    tok = _next_token(&cursor);
    if (!tok) {
        OSH_DIAG_ERRORF(diag,
                        "%s line %d: DCM card requires a patient position token "
                        "(e.g. HFS, HFP, FFS, FFP, HFDL, HFDR, FFDL, FFDR)",
                        geo_filename,
                        lineno);
        goto done;
    }
    patient_pos = osh_patient_position_from_str(tok);
    if (patient_pos == OSH_PP_UNKNOWN) {
        OSH_DIAG_ERRORF(diag,
                        "%s line %d: unknown DCM patient position '%s' "
                        "(expected one of: HFS HFP FFS FFP HFDL HFDR FFDL FFDR)",
                        geo_filename,
                        lineno,
                        tok);
        goto done;
    }

    tok = _next_token(&cursor);
    if (!_parse_double_token(tok, &gantry_deg)) {
        OSH_DIAG_ERRORF(diag, "%s line %d: invalid DCM gantry angle", geo_filename, lineno);
        goto done;
    }
    tok = _next_token(&cursor);
    if (!_parse_double_token(tok, &couch_deg)) {
        OSH_DIAG_ERRORF(diag, "%s line %d: invalid DCM couch angle", geo_filename, lineno);
        goto done;
    }
    tok = _next_token(&cursor);
    if (!_parse_double_token(tok, &iso_mm[0])) {
        OSH_DIAG_ERRORF(diag, "%s line %d: invalid DCM isocenter X (DICOM LPS mm)", geo_filename, lineno);
        goto done;
    }
    tok = _next_token(&cursor);
    if (!_parse_double_token(tok, &iso_mm[1])) {
        OSH_DIAG_ERRORF(diag, "%s line %d: invalid DCM isocenter Y (DICOM LPS mm)", geo_filename, lineno);
        goto done;
    }
    tok = _next_token(&cursor);
    if (!_parse_double_token(tok, &iso_mm[2])) {
        OSH_DIAG_ERRORF(diag, "%s line %d: invalid DCM isocenter Z (DICOM LPS mm)", geo_filename, lineno);
        goto done;
    }
    if (_next_token(&cursor) != NULL) {
        OSH_DIAG_ERRORF(diag,
                        "%s line %d: too many DCM arguments "
                        "(expected: name dir patient_pos gantry couch iso_x_mm iso_y_mm iso_z_mm)",
                        geo_filename,
                        lineno);
        goto done;
    }

    geo_dir = osh_path_dirname(geo_filename);
    if (osh_relative_path_to_file(&ct_dir_abs, geo_dir, ct_dir_rel) != 0) {
        rc = OSH_ENOMEM;
        goto done;
    }
    osh_path_normalize(ct_dir_abs);

    rc = osh_dicom_ct_read(ct_dir_abs, &ct, diag);
    if (rc != OSH_OK) {
        OSH_DIAG_ERRORF(diag, "%s line %d: failed to read CT DICOM directory '%s'", geo_filename, lineno, ct_dir_abs);
        goto done;
    }
    if (ct.rows <= 0 || ct.cols <= 0 || ct.n_slices <= 0 || ct.pixel_spacing[0] <= 0.0 || ct.pixel_spacing[1] <= 0.0
        || ct.slice_spacing <= 0.0 || !ct.pixels) {
        rc = OSH_EPARSE;
        OSH_DIAG_ERRORF(diag, "%s line %d: DCM CT metadata is invalid for voxel body setup", geo_filename, lineno);
        goto done;
    }
    if (!_ct_orientation_is_axial(&ct)) {
        rc = OSH_EPARSE;
        OSH_DIAG_ERRORF(diag,
                        "%s line %d: non-axial CT orientation is not supported yet in DCM placement "
                        "(future work: apply row/col cosines); row=[%.6f %.6f %.6f], col=[%.6f %.6f %.6f]",
                        geo_filename,
                        lineno,
                        ct.row_cosine[0],
                        ct.row_cosine[1],
                        ct.row_cosine[2],
                        ct.col_cosine[0],
                        ct.col_cosine[1],
                        ct.col_cosine[2]);
        goto done;
    }

    /* Compute base rotation for the chosen patient position.
     * tb[i] is the i-th DICOM axis expressed in universe coordinates.
     * The mapping is DICOM_local[i] = tb[i] . p_universe.
     * This embeds the IEC 61217 patient orientation into the transform chain. */
    osh_patient_position_base_rotation(patient_pos, tb);

    /* Local voxel grid in cm plus transform terms consumed by _setup_vox().
     * Keep x0/y0/z0 explicit even for DCM: they are the generic VOX corner
     * parameters and may be non-zero for non-DCM VOX sources. */
    par_out[0] = 0.0; /* x0: local corner of voxel [0,0,0] */
    par_out[1] = 0.0; /* y0: local corner of voxel [0,0,0] */
    par_out[2] = 0.0; /* z0: local corner of voxel [0,0,0] */
    par_out[3] = 0.1 * ct.pixel_spacing[1];
    par_out[4] = 0.1 * ct.pixel_spacing[0];
    par_out[5] = 0.1 * ct.slice_spacing;
    par_out[6] = (double) ct.cols;
    par_out[7] = (double) ct.rows;
    par_out[8] = (double) ct.n_slices;
    par_out[9] = gantry_deg;
    par_out[10] = couch_deg;

    /* Apply couch then gantry rotations to tb (same order as _setup_vox / _vox_body_build_transform). */
    gantry_rad = gantry_deg * OSH_M_PI_180;
    couch_rad = couch_deg * OSH_M_PI_180;
    for (j = 0; j < 3; j++) {
        osh_vect_rot_z(-couch_rad, tb[j]); /* IEC couch: CCW from above; rot_z is CW, so negate */
        osh_vect_rot_y(gantry_rad, tb[j]);
    }

    /* Compute isocenter in CT-local frame (cm, relative to CT corner).
     * ct.origin[j] is the DICOM LPS position of the first voxel CENTER (mm).
     * CT-local axis j corresponds to DICOM axis j with spacing par_out[3+j] cm.
     * Adding spacing/2 converts DICOM center-based origin to corner-based local origin. */
    for (j = 0; j < 3; j++) {
        iso_ct_local[j] = (iso_mm[j] - ct.origin[j]) * 0.1 + par_out[3 + j] * 0.5;
    }

    /* t_corner = -R^T * iso_ct_local, where R = rotated tb matrix.
     * The body transform maps universe -> local via p_local = R*(p_universe - t_corner),
     * so t_corner is the universe position that maps to CT local origin (0,0,0).
     * At gantry/couch=0 this equals the CT corner in universe; at non-zero angles
     * it is computed correctly for any orientation.
     * R^T[k][j] = tb[j][k], so t_corner[k] = -sum_j(tb[j][k] * iso_ct_local[j]). */
    for (k = 0; k < 3; k++) {
        par_out[11 + k] = 0.0;
        for (j = 0; j < 3; j++) {
            par_out[11 + k] -= tb[j][k] * iso_ct_local[j];
        }
    }
    /* CT DICOM origin offset [cm]: par_out[14+j] = -ct.origin[j]/10.
     *
     * WHY: RTDOSE/DicomCT scoring grids are expressed in CT-local coordinates,
     * i.e. DICOM physical space but relative to the CT's first voxel corner.
     * The RTDOSE DICOM origin is an absolute DICOM coordinate (mm); adding
     * par_out[14+j] converts it to CT-local:
     *   CT_local[j] = (rd.origin[j] - ct.origin[j]) / 10
     *               = rd.origin[j]/10 + par_out[14+j]
     *
     * No patient-position or gantry rotation terms appear here because the
     * RTDOSE grid is already expressed in DICOM physical (LPS) space, which is
     * the same coordinate system as the CT voxel grid.  The universe<->DICOM
     * rotation is handled entirely at scoring time via the has_rotation/ct_t
     * path in run_setup_voxel_scoring(), not at grid setup time. */
    for (j = 0; j < 3; j++) {
        par_out[14 + j] = -ct.origin[j] * 0.1;
    }
    /* Store patient position code for _setup_vox() and run_setup_voxel_scoring(). */
    par_out[17] = (double) (int) patient_pos;
    *npar_out = 18;
    if (OSH_VOXEL_LAYOUT_DEFAULT == OSH_VOXEL_ORDER_ROW_MAJOR) {
        size_t nxy;
        if (_mul_size_overflow((size_t) ct.cols, (size_t) ct.rows, &nxy)
            || _mul_size_overflow(nxy, (size_t) ct.n_slices, n_hu_out)) {
            rc = OSH_ENOMEM;
            goto done;
        }
        /* Row-major baseline: transfer ownership directly and avoid a full-volume copy. */
        *hu_out = ct.pixels;
        ct.pixels = NULL;
    } else {
        *hu_out = (int16_t *) osh_voxel_reorder(ct.pixels,
                                                (size_t) ct.cols,
                                                (size_t) ct.rows,
                                                (size_t) ct.n_slices,
                                                sizeof(int16_t),
                                                OSH_VOXEL_LAYOUT_DEFAULT,
                                                n_hu_out);
        free(ct.pixels);
        ct.pixels = NULL;
        if (!*hu_out) {
            rc = OSH_ENOMEM;
            goto done;
        }
    }
    rc = OSH_OK;

done:
    if (ct.pixels) {
        osh_dicom_ct_free(&ct);
    }
    free(ct_dir_abs);
    free(geo_dir);
    free(work);
    return rc;
}

/**
 * @brief Normalize legacy numeric material names.
 *
 * @details Maps "0" to "blackhole" and "1000" to "vacuum" with a warning;
 * all other names pass through unchanged.
 *
 * @param[in] raw       Raw material name token from the file.
 * @param[in] filename  Source file name used in warning messages.
 * @param[in] lineno    Source line number used in warning messages.
 *
 * @returns Normalized name string (either @p raw or a literal constant).
 */
static char const *
_normalize_material_name(char const *raw, char const *filename, int lineno, struct osh_diag_sink const *diag) {
    if (strcmp(raw, "0") == 0) {
        OSH_DIAG_WARNF(diag,
                       "%s line %d: legacy material '0' mapped to 'blackhole'; use 'blackhole' explicitly",
                       filename,
                       lineno);
        return "blackhole";
    }
    if (strcmp(raw, "1000") == 0) {
        OSH_DIAG_WARNF(
            diag, "%s line %d: legacy material '1000' mapped to 'vacuum'; use 'vacuum' explicitly", filename, lineno);
        return "vacuum";
    }
    return raw;
}

/**
 * @brief Duplicate a normalized material name into @p z->material_name.
 *
 * @param[in,out] z         Zone to update.
 * @param[in]     raw_name  Raw material name token (will be normalized).
 * @param[in]     filename  Source file name used in warning messages.
 * @param[in]     lineno    Source line number used in warning messages.
 *
 * @returns OSH_OK on success, OSH_ENOMEM on allocation failure.
 */
static enum osh_status _assign_material(struct osh_geometry_zone *z,
                                        char const *raw_name,
                                        char const *filename,
                                        int lineno,
                                        struct osh_diag_sink const *diag) {
    char const *name = _normalize_material_name(raw_name, filename, lineno, diag);
    char *copy = strdup(name);
    if (!copy) {
        return OSH_ENOMEM;
    }
    free(z->material_name);
    z->material_name = copy;
    return OSH_OK;
}

/**
 * @brief Rewind @p shf to the beginning and reset its line counter.
 *
 * @param[in,out] shf  Open geometry file to rewind.
 *
 * @returns 1 on success, 0 on I/O error.
 */
static int _rewind_oshfile(struct oshfile *shf, struct osh_diag_sink const *diag) {
    if (fseek(shf->fp, 0L, SEEK_SET) != 0) {
        OSH_DIAG_ERRORF(diag, "Failed to rewind geometry file '%s'", shf->filename);
        return 0;
    }
    shf->lineno = 0;
    return 1;
}

/* ---- Phase 1: body section ----------------------------------------------- */

/**
 * @brief Parse the body section from file start up to the first END card.
 *
 * @details Fills ws->bodies[0..ws->nbodies-1] with body type, name, and raw
 * argument arrays. Leaves the file positioned just after the first END card.
 *
 * @param[in,out] shf  Open geometry file; rewound internally.
 * @param[in,out] ws   Geometry workspace with pre-allocated bodies array.
 *
 * @returns OSH_OK on success, OSH_EPARSE or OSH_ENOMEM on failure.
 */
static enum osh_status
_parse_bodies(struct oshfile *shf, struct osh_diag_sink const *diag, struct osh_geometry_workspace *ws) {
    char *line = NULL;
    char *key = NULL;
    char *args = NULL;
    int lineno;
    int btype_new;
    int btype = OSH_GEOMETRY_BODY_NONE;
    int nt;
    int npar = 0;
    int off = 0;
    int allow_continuation = 0;
    char nstr[OSH_GEMCA_BODY_NAME_MAXLEN];
    double par[OSH_GEMCA_NARGS_MAX];
    int16_t *body_hu = NULL;
    size_t body_n_hu = 0u;
    size_t ibody = 0;
    int body_active = 0;
    enum osh_status rc = OSH_OK;

    if (!_rewind_oshfile(shf, diag)) {
        return OSH_EIO;
    }

    while (osh_readline_key(shf, &line, &key, &args, &lineno) > 0) {

        if (strcasecmp(key, OSH_GEO_KEY_END) == 0) {
            /* Finalize last body. */
            if (body_active) {
                rc = _finalize_body(ws, ibody, btype, nstr, par, npar);
                if (rc != OSH_OK) {
                    free(body_hu);
                    body_hu = NULL;
                    goto done;
                }
                ws->bodies[ibody].hu = body_hu;
                ws->bodies[ibody].n_hu = body_n_hu;
                body_hu = NULL;
            }
            free(line);
            line = NULL;
            return OSH_OK; /* File is now positioned just after the first END. */
        }

        btype_new = _body_type_from_key(key);

        if (btype_new != OSH_GEOMETRY_BODY_NONE) {
            /* Finalize the previous body, if any. */
            if (body_active) {
                rc = _finalize_body(ws, ibody, btype, nstr, par, npar);
                if (rc != OSH_OK) {
                    free(body_hu);
                    body_hu = NULL;
                    goto done;
                }
                ws->bodies[ibody].hu = body_hu;
                ws->bodies[ibody].n_hu = body_n_hu;
                body_hu = NULL;
                ibody++;
            }

            if (ibody >= ws->nbodies) {
                OSH_DIAG_ERRORF(diag, "%s line %d: too many bodies (max=%zu)", shf->filename, lineno, ws->nbodies);
                rc = OSH_EPARSE;
                goto done;
            }

            btype = btype_new;
            body_active = 1;

            if (!args) {
                OSH_DIAG_ERRORF(diag, "%s line %d: missing body name or parameters", shf->filename, lineno);
                rc = OSH_EPARSE;
                goto done;
            }
            allow_continuation = 1;

            if (strcasecmp(key, OSH_GEO_KEY_VOX) == 0) {
                OSH_DIAG_ERRORF(diag, "%s line %d: legacy VOX card is not supported; use DCM", shf->filename, lineno);
                rc = OSH_EPARSE;
                goto done;
            }
            if (strcasecmp(key, OSH_GEO_KEY_DCM) == 0) {
                int16_t *dcm_hu = NULL;
                size_t dcm_n_hu = 0u;
                rc = _parse_dcm_body(args, shf->filename, lineno, diag, nstr, par, &npar, &dcm_hu, &dcm_n_hu);
                if (rc != OSH_OK) {
                    free(dcm_hu);
                    goto done;
                }
                body_hu = dcm_hu;
                body_n_hu = dcm_n_hu;
                allow_continuation = 0;
                off = npar;
            } else {
                /* First six values: name + up to 5 floats on the same line. */
                nt = sscanf(
                    args, "%s %lf %lf %lf %lf %lf %lf", nstr, &par[0], &par[1], &par[2], &par[3], &par[4], &par[5]);
                if (nt < 1) {
                    OSH_DIAG_ERRORF(diag, "%s line %d: missing body name", shf->filename, lineno);
                    rc = OSH_EPARSE;
                    goto done;
                }
                npar = nt - 1;
                off = 6;
                body_hu = NULL;
                body_n_hu = 0u;
            }

        } else {
            /* Continuation line: accumulate more float arguments. */
            if (!body_active) {
                /* Allow preamble/header lines before the first body card. */
                free(line);
                line = NULL;
                continue;
            }
            if (!allow_continuation) {
                OSH_DIAG_ERRORF(
                    diag, "%s line %d: unexpected continuation line for this body card", shf->filename, lineno);
                rc = OSH_EPARSE;
                goto done;
            }

            if ((off + 5) >= OSH_GEMCA_NARGS_MAX) {
                OSH_DIAG_ERRORF(diag, "%s line %d: too many body arguments", shf->filename, lineno);
                rc = OSH_EPARSE;
                goto done;
            }

            nt = sscanf(key, "%lf", &par[off]);
            if (nt > 0) {
                npar += nt;
            }
            nt = 0;
            if (args) {
                nt = sscanf(args,
                            "%lf %lf %lf %lf %lf",
                            &par[1 + off],
                            &par[2 + off],
                            &par[3 + off],
                            &par[4 + off],
                            &par[5 + off]);
            }
            if (nt > 0) {
                npar += nt;
            }
            off += 6;
        }

        free(line);
        line = NULL;
    }

done:
    free(body_hu);
    free(line);
    return rc;
}

/* ---- Phase 2: zone section ----------------------------------------------- */

/**
 * @brief Parse the zone section between the first and second END cards.
 *
 * @details Fills ws->zones[].name and ws->zones[].expr. Precondition: file is
 * positioned just after the first END card as left by _parse_bodies().
 *
 * @param[in,out] shf  Open geometry file; read from current position.
 * @param[in,out] ws   Geometry workspace with pre-allocated zones array.
 *
 * @returns OSH_OK on success, OSH_EPARSE or OSH_ENOMEM on failure.
 */
static enum osh_status
_parse_zones(struct oshfile *shf, struct osh_diag_sink const *diag, struct osh_geometry_workspace *ws) {
    char *line = NULL;
    char *key = NULL;
    char *args = NULL;
    int lineno;
    size_t izone = 0u;
    int zone_active = 0;
    char *expr = NULL;
    enum osh_status rc = OSH_OK;

    /* Accumulator for the current zone's raw boolean expression.
     * Starts as an empty NUL-terminated string; grows via _str_append(). */
    expr = (char *) calloc(1u, sizeof(char));
    if (!expr) {
        return OSH_ENOMEM;
    }

    while (osh_readline_key(shf, &line, &key, &args, &lineno) > 0) {

        if (_is_zone_continuation(key)) {
            if (!zone_active) {
                OSH_DIAG_ERRORF(diag, "%s line %d: zone continuation before any zone header", shf->filename, lineno);
                rc = OSH_EPARSE;
                goto done;
            }
            rc = _str_append(&expr, key);
            if (rc != OSH_OK) {
                goto done;
            }

        } else {
            /* Non-continuation key: either a new zone name or the END card.
             * In both cases, finalize the expression for the current zone first. */
            if (zone_active && expr[0] != '\0') {
                if (izone >= ws->nzones) {
                    OSH_DIAG_ERRORF(diag, "%s line %d: too many zones (max=%zu)", shf->filename, lineno, ws->nzones);
                    rc = OSH_EPARSE;
                    goto done;
                }
                ws->zones[izone].expr = strdup(expr);
                if (!ws->zones[izone].expr) {
                    rc = OSH_ENOMEM;
                    goto done;
                }
                expr[0] = '\0'; /* reset accumulator for next zone */
            }

            if (strcasecmp(key, OSH_GEO_KEY_END) == 0) {
                free(line);
                line = NULL;
                goto done; /* File now positioned just after the second END. */
            }

            /* New zone: advance the index (except for the very first zone). */
            if (!zone_active) {
                zone_active = 1;
                /* izone stays at 0 */
            } else {
                izone++;
            }

            if (izone >= ws->nzones) {
                OSH_DIAG_ERRORF(diag, "%s line %d: too many zones (max=%zu)", shf->filename, lineno, ws->nzones);
                rc = OSH_EPARSE;
                goto done;
            }

            ws->zones[izone].name = strdup(key);
            if (!ws->zones[izone].name) {
                rc = OSH_ENOMEM;
                goto done;
            }
        }

        /* Append args (if any) to the current zone's expression accumulator. */
        if (args) {
            rc = _str_append(&expr, args);
            if (rc != OSH_OK) {
                goto done;
            }
        }

        free(line);
        line = NULL;
    }

done:
    free(line);
    free(expr);
    return rc;
}

/* ---- Phase 3: material section ------------------------------------------- */

/**
 * @brief Handle one ASSIGNMAT / ASSIGNMA card.
 *
 * @details Assigns a material to a range of zones specified by name, with an
 * optional stride. Updates ws->zones[].material_name for each matching zone.
 *
 * @param[in,out] ws        Geometry workspace whose zones are updated.
 * @param[in,out] args      Argument string after the ASSIGNMAT keyword; tokenized in place.
 * @param[in]     filename  Source file name used in error messages.
 * @param[in]     lineno    Source line number used in error messages.
 *
 * @returns OSH_OK on success, OSH_EPARSE on invalid arguments, OSH_ENOMEM on failure.
 */
static enum osh_status _parse_assignmat(
    struct osh_geometry_workspace *ws, char *args, char const *filename, int lineno, struct osh_diag_sink const *diag) {
    char *cursor = args;
    char *mat_name;
    char *zname_start;
    char *zname_end;
    char *stride_str;
    size_t iz_start = 0u;
    size_t iz_end = 0u;
    size_t stride;
    size_t iz;
    enum osh_status rc;

    mat_name = _next_token(&cursor);
    if (!mat_name) {
        OSH_DIAG_ERRORF(diag, "%s line %d: ASSIGNMAT: missing material name", filename, lineno);
        return OSH_EPARSE;
    }

    zname_start = _next_token(&cursor);
    if (!zname_start) {
        OSH_DIAG_ERRORF(diag, "%s line %d: ASSIGNMAT: missing zone name", filename, lineno);
        return OSH_EPARSE;
    }
    if (!_zone_index_from_name(zname_start, ws->zones, ws->nzones, &iz_start)) {
        OSH_DIAG_ERRORF(diag, "%s line %d: ASSIGNMAT: unknown zone '%s'", filename, lineno, zname_start);
        return OSH_EPARSE;
    }

    zname_end = _next_token(&cursor);
    if (zname_end) {
        if (!_zone_index_from_name(zname_end, ws->zones, ws->nzones, &iz_end)) {
            OSH_DIAG_ERRORF(diag, "%s line %d: ASSIGNMAT: unknown zone '%s'", filename, lineno, zname_end);
            return OSH_EPARSE;
        }
    } else {
        iz_end = iz_start;
    }

    stride_str = _next_token(&cursor);
    stride = stride_str ? (size_t) atoi(stride_str) : 1u;
    if (stride == 0u) {
        OSH_DIAG_ERRORF(diag, "%s line %d: ASSIGNMAT: invalid stride 0", filename, lineno);
        return OSH_EPARSE;
    }

    if (iz_end < iz_start) {
        OSH_DIAG_ERRORF(diag, "%s line %d: ASSIGNMAT: zone range ends before it starts", filename, lineno);
        return OSH_EPARSE;
    }

    for (iz = iz_start; iz <= iz_end; iz += stride) {
        rc = _assign_material(&ws->zones[iz], mat_name, filename, lineno, diag);
        if (rc != OSH_OK) {
            return rc;
        }
    }
    return OSH_OK;
}

/**
 * @brief Parse the material section after the second END card.
 *
 * @details Handles both ASSIGNMAT cards and the legacy positional zone-list
 * format. Rewinds the file internally and skips past both END cards.
 *
 * @param[in,out] shf  Open geometry file; rewound internally.
 * @param[in,out] ws   Geometry workspace whose zones receive material names.
 *
 * @returns OSH_OK on success, OSH_EPARSE or OSH_ENOMEM on failure.
 */
static enum osh_status
_parse_media(struct oshfile *shf, struct osh_diag_sink const *diag, struct osh_geometry_workspace *ws) {
    char *line = NULL;
    char *key = NULL;
    char *args = NULL;
    int lineno;
    int nend = 0;
    int in_media = 0;
    int warned_legacy = 0;
    size_t izone = 0u;
    char *cursor;
    char *tok;
    enum osh_status rc = OSH_OK;

    if (!_rewind_oshfile(shf, diag)) {
        return OSH_EIO;
    }

    /* Skip past the second END card. */
    while (osh_readline_key(shf, &line, &key, &args, &lineno) > 0) {
        if (strcasecmp(key, OSH_GEO_KEY_END) == 0) {
            ++nend;
            if (nend == 2) {
                free(line);
                line = NULL;
                break;
            }
        }
        free(line);
        line = NULL;
    }

    if (nend < 2) {
        /* File has no material section — nothing to parse. */
        free(line);
        return OSH_OK;
    }

    /* Parse material assignments.
     *
     * Legacy positional format:
     *   - First nzones tokens form the zone-number list (ignored; zones are
     *     already named).
     *   - Next nzones tokens are material names assigned positionally.
     *
     * Preferred ASSIGNMAT format:
     *   ASSIGNMAT <material> <zone-start> [<zone-end> [<stride>]]
     */
    while (osh_readline_key(shf, &line, &key, &args, &lineno) > 0) {

        if ((strcasecmp(key, OSH_GEO_KEY_ASSIGNMAT) == 0) || (strcasecmp(key, OSH_GEO_KEY_ASSIGNMA) == 0)) {
            rc = _parse_assignmat(ws, args, shf->filename, lineno, diag);
            if (rc != OSH_OK) {
                free(line);
                return rc;
            }
            free(line);
            line = NULL;
            continue;
        }

        /* Legacy: count tokens and assign materials when in the media half. */
        izone++;
        if (izone > ws->nzones) {
            OSH_DIAG_ERRORF(
                diag, "%s line %d: too many entries in material section (max=%zu)", shf->filename, lineno, ws->nzones);
            free(line);
            return OSH_EPARSE;
        }

        if (in_media) {
            if (!warned_legacy) {
                OSH_DIAG_WARNF(diag, "Implicit geo.dat material lists are legacy; use ASSIGNMAT instead");
                warned_legacy = 1;
            }
            rc = _assign_material(&ws->zones[izone - 1u], key, shf->filename, lineno, diag);
            if (rc != OSH_OK) {
                free(line);
                return rc;
            }
        }

        cursor = args;
        tok = _next_token(&cursor);
        while (tok != NULL) {
            izone++;
            if (izone > ws->nzones) {
                OSH_DIAG_ERRORF(diag,
                                "%s line %d: too many entries in material section (max=%zu)",
                                shf->filename,
                                lineno,
                                ws->nzones);
                free(line);
                return OSH_EPARSE;
            }
            if (in_media) {
                rc = _assign_material(&ws->zones[izone - 1u], tok, shf->filename, lineno, diag);
                if (rc != OSH_OK) {
                    free(line);
                    return rc;
                }
            }
            tok = _next_token(&cursor);
        }

        /* After counting nzones entries we switch from the zone-number list
         * to the material-name list. */
        if (izone == ws->nzones) {
            in_media = 1;
            izone = 0u;
        }

        free(line);
        line = NULL;
    }

    free(line);
    return OSH_OK;
}

/* ---- Public entry point -------------------------------------------------- */

enum osh_status
osh_geometry_parse(struct oshfile *oshf, struct osh_diag_sink const *diag, struct osh_geometry_workspace *ws) {
    size_t nbodies;
    size_t nzones;
    enum osh_status rc;

    if (!oshf || !ws) {
        return OSH_EINVAL;
    }

    /* Count bodies and zones without retaining any state. */
    nbodies = _count_bodies(oshf, diag);
    nzones = _count_zones(oshf, diag);

    if (nbodies == 0u || nzones == 0u) {
        OSH_DIAG_ERRORF(diag, "geometry: '%s' has no bodies or no zones", oshf->filename);
        return OSH_EPARSE;
    }

    ws->nbodies = nbodies;
    ws->bodies = (struct osh_geometry_body *) calloc(nbodies, sizeof(struct osh_geometry_body));
    if (!ws->bodies) {
        return OSH_ENOMEM;
    }

    ws->nzones = nzones;
    ws->zones = (struct osh_geometry_zone *) calloc(nzones, sizeof(struct osh_geometry_zone));
    if (!ws->zones) {
        return OSH_ENOMEM;
    }

    /* Phase 1: parse bodies.
     * Rewinds internally; leaves file positioned after the first END card. */
    rc = _parse_bodies(oshf, diag, ws);
    if (rc != OSH_OK) {
        return rc;
    }

    /* Phase 2: parse zones.
     * Continues from where _parse_bodies stopped;
     * leaves file positioned after the second END card. */
    rc = _parse_zones(oshf, diag, ws);
    if (rc != OSH_OK) {
        return rc;
    }

    /* Phase 3: parse material assignments.
     * Rewinds internally; skips past both END cards. */
    rc = _parse_media(oshf, diag, ws);
    if (rc != OSH_OK) {
        return rc;
    }

    return OSH_OK;
}
