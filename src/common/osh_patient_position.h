#ifndef OSH_PATIENT_POSITION_H
#define OSH_PATIENT_POSITION_H

/*
 * IEC 61217 patient position support for the OpenShieldHIT simulation universe.
 *
 * WHY THIS FILE EXISTS
 * ====================
 * The simulation universe axes follow IEC 61217 conventions:
 *   Universe X = patient left (IEC Xf)
 *   Universe Y = along table, cranial direction for head-first patients (IEC body axis)
 *   Universe Z = vertically above the patient = patient anterior for HFS supine (IEC Yf)
 *               = nozzle direction at gantry angle 0 degrees
 *
 * DICOM stores image data in LPS (Left-Posterior-Superior) patient coordinates:
 *   DICOM X = patient left    (L)
 *   DICOM Y = patient posterior (P)
 *   DICOM Z = patient superior / cranial (S)
 *
 * The patient position (how the patient lies on the table) determines the
 * relationship between the DICOM patient axes and the IEC universe axes.
 * For example, a Head-First Supine (HFS) patient has:
 *   - Patient left = universe X  (same)
 *   - Patient cranial = universe Y, but cranial = +Z in DICOM -> so DICOM Z -> universe Y
 *   - Patient anterior (belly-up) = universe Z, but anterior = -Y in DICOM -> DICOM -Y -> universe Z
 *
 * The matrix tb[3][3] encodes this mapping: each row i gives the DICOM axis i
 * expressed as a unit vector in universe coordinates.  Concretely:
 *
 *   DICOM_local[i] = tb[i] . p_universe      (dot product)
 *
 * which means tb maps universe coordinates to DICOM patient coordinates.
 * The matrix is orthonormal (det=+1) for all standard positions.
 *
 * COORDINATE CHAIN SUMMARY
 * =========================
 *  1. Start with universe position p_uni (IEC 61217 frame).
 *  2. Apply tb (patient-position base rotation): DICOM_local[i] = tb[i] . p_uni
 *  3. Apply gantry rotation around universe Y (sagittal plane, IEC convention).
 *  4. Apply couch  rotation around universe Z (vertical axis, IEC convention).
 *  The combined transform maps universe positions to CT voxel-local positions.
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief IEC 61217 patient position codes.
 *
 * @details
 * Standard patient positions as defined in IEC 61217 and used in DICOM
 * PatientPosition (0018,5100).  The integer values are stored in VOX body
 * parameter a[17] so they must remain stable across versions.
 *
 * OSH_PP_UNKNOWN (-1) is the sentinel for an unrecognised position string.
 * OSH_PP_HFS (0) is used as the default when a[17] is absent (old geo.dat
 * files without the patient-position token).
 */
enum osh_patient_position {
    OSH_PP_HFS = 0,  /* Head First Supine       — standard brain/thorax position */
    OSH_PP_HFP = 1,  /* Head First Prone        */
    OSH_PP_FFS = 2,  /* Feet First Supine       */
    OSH_PP_FFP = 3,  /* Feet First Prone        */
    OSH_PP_HFDL = 4, /* Head First Decubitus Left  (lying on left side) */
    OSH_PP_HFDR = 5, /* Head First Decubitus Right (lying on right side) */
    OSH_PP_FFDL = 6, /* Feet First Decubitus Left  */
    OSH_PP_FFDR = 7, /* Feet First Decubitus Right */
    OSH_PP_UNKNOWN = -1
};

/**
 * @brief Parse a DICOM PatientPosition string to the corresponding enum value.
 *
 * @details
 * Accepts the standard DICOM PatientPosition (0018,5100) strings: "HFS",
 * "HFP", "FFS", "FFP", "HFDL", "HFDR", "FFDL", "FFDR".  Matching is
 * case-insensitive: the input is lowercased into a local buffer before
 * comparison with lowercase literals.
 *
 * WHY NOT strcasecmp: Windows (MSVC) does not provide strcasecmp.  The
 * project policy (see DEVELOPER.md) requires lowercasing at match time with a
 * manual tolower loop, then plain strcmp.
 *
 * WHY A SEPARATE .c FILE: keeping this function in a .c file (rather than
 * static-in-header) avoids -Wunused-function warnings in translation units
 * that include this header but only need osh_patient_position_base_rotation.
 *
 * @param[in] s  Patient position string, e.g. "HFS" or "hfs".
 *
 * @returns The matching enum value, or OSH_PP_UNKNOWN if @p s is NULL,
 *          empty, or not a recognised position string.
 */
enum osh_patient_position osh_patient_position_from_str(char const *s);

/**
 * @brief Fill the base rotation matrix for a given IEC 61217 patient position.
 *
 * @details
 * Fills @p tb[3][3] with the base rotation matrix that maps universe
 * coordinates to DICOM LPS patient coordinates.  Each row i encodes:
 *
 *   DICOM_axis_i_in_universe = tb[i]   (the i-th DICOM axis as a universe vector)
 *
 * Equivalently: DICOM_local[i] = tb[i] . p_universe.
 *
 * IEC 61217 UNIVERSE AXES (reminder):
 *   Universe X = patient left
 *   Universe Y = cranial direction for head-first patients
 *   Universe Z = patient anterior for HFS supine = nozzle direction at gantry 0
 *
 * DICOM LPS AXES (reminder):
 *   DICOM X = patient left     (L)
 *   DICOM Y = patient posterior (P)
 *   DICOM Z = patient superior  (S, i.e. cranial)
 *
 * The rows of tb express each DICOM axis as a vector in universe coordinates:
 *
 *   tb[0] = DICOM X in universe = patient-left direction
 *   tb[1] = DICOM Y in universe = patient-posterior direction
 *   tb[2] = DICOM Z in universe = patient-cranial direction
 *
 * DERIVATION FOR HFS (Head First Supine) — the most common case:
 *   Patient lies on back, head toward gantry (cranial = +Y universe).
 *   Patient left  = universe X  -> tb[0] = [1, 0, 0]
 *   Patient post. = universe -Z (supine: back is down, -Z = below patient)
 *                                -> tb[1] = [0, 0, -1]
 *     (Universe Z points anteriorly/upward from patient; posterior = opposite)
 *   Patient cranial = universe Y -> tb[2] = [0, 1, 0]
 *
 * DERIVATION FOR HFP (Head First Prone):
 *   Patient lies face-down, head toward gantry.
 *   Patient left  becomes -universe X (patient is flipped L-R when prone)
 *                                -> tb[0] = [-1, 0, 0]
 *   Patient post. = universe +Z (belly faces down, posterior = up = +Z)
 *                                -> tb[1] = [0, 0, 1]
 *   Patient cranial = universe Y -> tb[2] = [0, 1, 0]
 *
 * DERIVATION FOR FFS (Feet First Supine):
 *   Patient lies on back, feet toward gantry (cranial = -Y universe).
 *   Patient left  = -universe X (feet-first flips L-R)
 *                                -> tb[0] = [-1, 0, 0]
 *   Patient post. = universe -Z (supine: back down) -> tb[1] = [0, 0, -1]
 *   Patient cranial = -universe Y -> tb[2] = [0, -1, 0]
 *
 * DERIVATION FOR FFP (Feet First Prone):
 *   Patient lies face-down, feet toward gantry.
 *   Patient left  = universe X  -> tb[0] = [1, 0, 0]
 *   Patient post. = universe +Z -> tb[1] = [0, 0, 1]
 *   Patient cranial = -universe Y -> tb[2] = [0, -1, 0]
 *
 * DECUBITUS POSITIONS: patient lies on their side.
 *   HFDL (Head First Decubitus Left — patient's LEFT side is down):
 *     Patient left  = -universe Z (left side faces down toward table = -Z)
 *                                  -> tb[0] = [0, 0, -1]
 *     Patient post. = -universe X (posterior faces opposite to patient left)
 *                                  -> tb[1] = [-1, 0, 0]
 *     Patient cranial = universe Y -> tb[2] = [0, 1, 0]
 *
 *   HFDR (Head First Decubitus Right — patient's RIGHT side is down):
 *     Patient left  = +universe Z -> tb[0] = [0, 0, 1]
 *     Patient post. = +universe X -> tb[1] = [1, 0, 0]
 *     Patient cranial = universe Y -> tb[2] = [0, 1, 0]
 *
 *   FFDL (Feet First Decubitus Left):
 *     Patient left  = -universe Z -> tb[0] = [0, 0, -1]
 *     Patient post. = +universe X -> tb[1] = [1, 0, 0]
 *     Patient cranial = -universe Y -> tb[2] = [0, -1, 0]
 *
 *   FFDR (Feet First Decubitus Right):
 *     Patient left  = +universe Z -> tb[0] = [0, 0, 1]
 *     Patient post. = -universe X -> tb[1] = [-1, 0, 0]
 *     Patient cranial = -universe Y -> tb[2] = [0, -1, 0]
 *
 * All eight matrices have determinant +1 (proper rotation, no reflection).
 * HFS is the default when no patient position is specified (legacy behaviour).
 *
 * @param[in]  pp  Patient position enum value.
 * @param[out] tb  3x3 matrix filled with the base rotation (row i = DICOM
 *                 axis i expressed in universe coordinates).
 */
void osh_patient_position_base_rotation(enum osh_patient_position pp, double tb[3][3]);

#ifdef __cplusplus
}
#endif

#endif /* OSH_PATIENT_POSITION_H */
