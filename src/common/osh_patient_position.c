#include "common/osh_patient_position.h"

#include <string.h>

/**
 * @brief Fill the base rotation matrix for a given IEC 61217 patient position.
 *
 * @details
 * See header for full documentation, derivations, and sign conventions.
 *
 * Each row of tb encodes a DICOM LPS axis expressed in universe coordinates:
 *   tb[i] = DICOM axis i in universe coords
 *   DICOM_local[i] = tb[i] . p_universe
 *
 * All eight matrices have determinant +1 (proper rotation, no reflection).
 * The default (unknown/fallback) is HFS, which is the most common clinical
 * position and keeps legacy geo.dat files (without patient_pos token) valid.
 */
void osh_patient_position_base_rotation(enum osh_patient_position pp, double tb[3][3]) {
    int i;
    int j;

    /* Zero-initialise so only the non-zero elements below need to be set. */
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            tb[i][j] = 0.0;
        }
    }

    switch (pp) {
        case OSH_PP_HFS:
            /* Head First Supine: the most common clinical position.
             * DICOM X = patient left  = universe  X
             * DICOM Y = patient post. = universe -Z  (supine: back is down, -Z = toward table)
             * DICOM Z = patient cran. = universe  Y  (head toward gantry = +Y) */
            tb[0][0] =  1.0;  /* DICOM X -> universe +X */
            tb[1][2] = -1.0;  /* DICOM Y -> universe -Z */
            tb[2][1] =  1.0;  /* DICOM Z -> universe +Y */
            break;

        case OSH_PP_HFP:
            /* Head First Prone: patient face-down, head toward gantry.
             * L-R is mirrored relative to HFS (left/right swap when prone).
             * DICOM X = patient left  = universe -X
             * DICOM Y = patient post. = universe +Z  (prone: back faces up = +Z)
             * DICOM Z = patient cran. = universe  Y */
            tb[0][0] = -1.0;  /* DICOM X -> universe -X */
            tb[1][2] =  1.0;  /* DICOM Y -> universe +Z */
            tb[2][1] =  1.0;  /* DICOM Z -> universe +Y */
            break;

        case OSH_PP_FFS:
            /* Feet First Supine: patient on back, feet toward gantry.
             * Feet-first swaps both L-R and cranial/caudal relative to HFS.
             * DICOM X = patient left  = universe -X (L-R flip for feet-first)
             * DICOM Y = patient post. = universe -Z  (supine: back down)
             * DICOM Z = patient cran. = universe -Y  (cranial now points in -Y) */
            tb[0][0] = -1.0;  /* DICOM X -> universe -X */
            tb[1][2] = -1.0;  /* DICOM Y -> universe -Z */
            tb[2][1] = -1.0;  /* DICOM Z -> universe -Y */
            break;

        case OSH_PP_FFP:
            /* Feet First Prone: patient face-down, feet toward gantry.
             * DICOM X = patient left  = universe +X
             * DICOM Y = patient post. = universe +Z  (prone: back up)
             * DICOM Z = patient cran. = universe -Y */
            tb[0][0] =  1.0;  /* DICOM X -> universe +X */
            tb[1][2] =  1.0;  /* DICOM Y -> universe +Z */
            tb[2][1] = -1.0;  /* DICOM Z -> universe -Y */
            break;

        case OSH_PP_HFDL:
            /* Head First Decubitus Left: patient's LEFT side faces down toward table.
             * DICOM X = patient left  = universe -Z  (left side pressed to table = -Z)
             * DICOM Y = patient post. = universe -X
             * DICOM Z = patient cran. = universe  Y */
            tb[0][2] = -1.0;  /* DICOM X -> universe -Z */
            tb[1][0] = -1.0;  /* DICOM Y -> universe -X */
            tb[2][1] =  1.0;  /* DICOM Z -> universe +Y */
            break;

        case OSH_PP_HFDR:
            /* Head First Decubitus Right: patient's RIGHT side faces down toward table.
             * DICOM X = patient left  = universe +Z
             * DICOM Y = patient post. = universe +X
             * DICOM Z = patient cran. = universe  Y */
            tb[0][2] =  1.0;  /* DICOM X -> universe +Z */
            tb[1][0] =  1.0;  /* DICOM Y -> universe +X */
            tb[2][1] =  1.0;  /* DICOM Z -> universe +Y */
            break;

        case OSH_PP_FFDL:
            /* Feet First Decubitus Left: patient's LEFT side down, feet toward gantry.
             * DICOM X = patient left  = universe -Z
             * DICOM Y = patient post. = universe +X
             * DICOM Z = patient cran. = universe -Y */
            tb[0][2] = -1.0;  /* DICOM X -> universe -Z */
            tb[1][0] =  1.0;  /* DICOM Y -> universe +X */
            tb[2][1] = -1.0;  /* DICOM Z -> universe -Y */
            break;

        case OSH_PP_FFDR:
            /* Feet First Decubitus Right: patient's RIGHT side down, feet toward gantry.
             * DICOM X = patient left  = universe +Z
             * DICOM Y = patient post. = universe -X
             * DICOM Z = patient cran. = universe -Y */
            tb[0][2] =  1.0;  /* DICOM X -> universe +Z */
            tb[1][0] = -1.0;  /* DICOM Y -> universe -X */
            tb[2][1] = -1.0;  /* DICOM Z -> universe -Y */
            break;

        default:
            /* Unrecognised code: fall back to HFS to keep transport valid.
             * The caller should have rejected OSH_PP_UNKNOWN at parse time;
             * this branch guards against accidental integer cast errors. */
            tb[0][0] =  1.0;
            tb[1][2] = -1.0;
            tb[2][1] =  1.0;
            break;
    }
}

/**
 * @brief Parse a DICOM PatientPosition string to the corresponding enum value.
 *
 * @details
 * See header for full documentation and design rationale.
 * Implementation: lowercase the input into a local buffer, then use strcmp
 * with lowercase literals.  strcasecmp is banned (not available on MSVC).
 */
enum osh_patient_position osh_patient_position_from_str(char const *s) {
    /* Local buffer sized for the longest token "HFDR" (4 chars) + NUL. */
    char buf[8];
    size_t i;
    size_t len;

    if (!s) {
        return OSH_PP_UNKNOWN;
    }
    len = strlen(s);
    if (len == 0u || len >= sizeof(buf)) {
        return OSH_PP_UNKNOWN;
    }
    /* Lowercase into buf so we can use strcmp with lowercase literals.
     * This avoids strcasecmp which is banned on Windows (see DEVELOPER.md). */
    for (i = 0u; i < len; i++) {
        char c = s[i];
        if (c >= 'A' && c <= 'Z') {
            c = (char)(c + ('a' - 'A'));
        }
        buf[i] = c;
    }
    buf[len] = '\0';

    if (strcmp(buf, "hfs")  == 0) { return OSH_PP_HFS;  }
    if (strcmp(buf, "hfp")  == 0) { return OSH_PP_HFP;  }
    if (strcmp(buf, "ffs")  == 0) { return OSH_PP_FFS;  }
    if (strcmp(buf, "ffp")  == 0) { return OSH_PP_FFP;  }
    if (strcmp(buf, "hfdl") == 0) { return OSH_PP_HFDL; }
    if (strcmp(buf, "hfdr") == 0) { return OSH_PP_HFDR; }
    if (strcmp(buf, "ffdl") == 0) { return OSH_PP_FFDL; }
    if (strcmp(buf, "ffdr") == 0) { return OSH_PP_FFDR; }
    return OSH_PP_UNKNOWN;
}
