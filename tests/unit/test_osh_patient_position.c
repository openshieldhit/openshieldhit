#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "common/osh_patient_position.h"
#include "test_assert.h"

static void test_from_str_all_positions(void);
static void test_from_str_case_insensitive(void);
static void test_from_str_invalid(void);
static void test_base_rotation_values(void);
static void test_base_rotation_orthonormal(void);
static void test_base_rotation_default_fallback(void);

static int nearly(double a, double b) {
    return fabs(a - b) < 1e-12;
}

/* Check that a 3x3 matrix is a proper rotation (orthonormal, det=+1). */
static void assert_proper_rotation(double tb[3][3]) {
    int i;
    int j;
    for (i = 0; i < 3; i++) {
        double norm2 = 0.0;
        for (j = 0; j < 3; j++)
            norm2 += tb[i][j] * tb[i][j];
        ASSERT_TRUE(nearly(norm2, 1.0));
    }
    for (i = 0; i < 3; i++) {
        for (j = i + 1; j < 3; j++) {
            double dot = tb[i][0] * tb[j][0] + tb[i][1] * tb[j][1] + tb[i][2] * tb[j][2];
            ASSERT_TRUE(nearly(dot, 0.0));
        }
    }
    /* det = +1 */
    double det = tb[0][0] * (tb[1][1] * tb[2][2] - tb[1][2] * tb[2][1])
                 - tb[0][1] * (tb[1][0] * tb[2][2] - tb[1][2] * tb[2][0])
                 + tb[0][2] * (tb[1][0] * tb[2][1] - tb[1][1] * tb[2][0]);
    ASSERT_TRUE(nearly(det, 1.0));
}

int main(void) {
    test_from_str_all_positions();
    test_from_str_case_insensitive();
    test_from_str_invalid();
    test_base_rotation_values();
    test_base_rotation_orthonormal();
    test_base_rotation_default_fallback();
    return 0;
}

static void test_from_str_all_positions(void) {
    ASSERT_TRUE(osh_patient_position_from_str("HFS") == OSH_PP_HFS);
    ASSERT_TRUE(osh_patient_position_from_str("HFP") == OSH_PP_HFP);
    ASSERT_TRUE(osh_patient_position_from_str("FFS") == OSH_PP_FFS);
    ASSERT_TRUE(osh_patient_position_from_str("FFP") == OSH_PP_FFP);
    ASSERT_TRUE(osh_patient_position_from_str("HFDL") == OSH_PP_HFDL);
    ASSERT_TRUE(osh_patient_position_from_str("HFDR") == OSH_PP_HFDR);
    ASSERT_TRUE(osh_patient_position_from_str("FFDL") == OSH_PP_FFDL);
    ASSERT_TRUE(osh_patient_position_from_str("FFDR") == OSH_PP_FFDR);
}

static void test_from_str_case_insensitive(void) {
    ASSERT_TRUE(osh_patient_position_from_str("hfs") == OSH_PP_HFS);
    ASSERT_TRUE(osh_patient_position_from_str("hfp") == OSH_PP_HFP);
    ASSERT_TRUE(osh_patient_position_from_str("ffs") == OSH_PP_FFS);
    ASSERT_TRUE(osh_patient_position_from_str("ffp") == OSH_PP_FFP);
    ASSERT_TRUE(osh_patient_position_from_str("hfdl") == OSH_PP_HFDL);
    ASSERT_TRUE(osh_patient_position_from_str("hfdr") == OSH_PP_HFDR);
    ASSERT_TRUE(osh_patient_position_from_str("ffdl") == OSH_PP_FFDL);
    ASSERT_TRUE(osh_patient_position_from_str("ffdr") == OSH_PP_FFDR);
    /* mixed case */
    ASSERT_TRUE(osh_patient_position_from_str("Hfs") == OSH_PP_HFS);
    ASSERT_TRUE(osh_patient_position_from_str("HfDr") == OSH_PP_HFDR);
}

static void test_from_str_invalid(void) {
    ASSERT_TRUE(osh_patient_position_from_str(NULL) == OSH_PP_UNKNOWN);
    ASSERT_TRUE(osh_patient_position_from_str("") == OSH_PP_UNKNOWN);
    ASSERT_TRUE(osh_patient_position_from_str("UNKNOWN") == OSH_PP_UNKNOWN);
    ASSERT_TRUE(osh_patient_position_from_str("HFSEXTRA") == OSH_PP_UNKNOWN); /* too long */
    ASSERT_TRUE(osh_patient_position_from_str("HF") == OSH_PP_UNKNOWN);
    ASSERT_TRUE(osh_patient_position_from_str("supine") == OSH_PP_UNKNOWN);
}

static void test_base_rotation_values(void) {
    double tb[3][3];

    /* HFS: DICOM X=+uX, DICOM Y=-uZ, DICOM Z=+uY */
    osh_patient_position_base_rotation(OSH_PP_HFS, tb);
    ASSERT_TRUE(nearly(tb[0][0], 1.0) && nearly(tb[0][1], 0.0) && nearly(tb[0][2], 0.0));
    ASSERT_TRUE(nearly(tb[1][0], 0.0) && nearly(tb[1][1], 0.0) && nearly(tb[1][2], -1.0));
    ASSERT_TRUE(nearly(tb[2][0], 0.0) && nearly(tb[2][1], 1.0) && nearly(tb[2][2], 0.0));

    /* HFP: DICOM X=-uX, DICOM Y=+uZ, DICOM Z=+uY */
    osh_patient_position_base_rotation(OSH_PP_HFP, tb);
    ASSERT_TRUE(nearly(tb[0][0], -1.0) && nearly(tb[0][1], 0.0) && nearly(tb[0][2], 0.0));
    ASSERT_TRUE(nearly(tb[1][0], 0.0) && nearly(tb[1][1], 0.0) && nearly(tb[1][2], 1.0));
    ASSERT_TRUE(nearly(tb[2][0], 0.0) && nearly(tb[2][1], 1.0) && nearly(tb[2][2], 0.0));

    /* FFS: DICOM X=-uX, DICOM Y=-uZ, DICOM Z=-uY */
    osh_patient_position_base_rotation(OSH_PP_FFS, tb);
    ASSERT_TRUE(nearly(tb[0][0], -1.0) && nearly(tb[0][1], 0.0) && nearly(tb[0][2], 0.0));
    ASSERT_TRUE(nearly(tb[1][0], 0.0) && nearly(tb[1][1], 0.0) && nearly(tb[1][2], -1.0));
    ASSERT_TRUE(nearly(tb[2][0], 0.0) && nearly(tb[2][1], -1.0) && nearly(tb[2][2], 0.0));

    /* FFP: DICOM X=+uX, DICOM Y=+uZ, DICOM Z=-uY */
    osh_patient_position_base_rotation(OSH_PP_FFP, tb);
    ASSERT_TRUE(nearly(tb[0][0], 1.0) && nearly(tb[0][1], 0.0) && nearly(tb[0][2], 0.0));
    ASSERT_TRUE(nearly(tb[1][0], 0.0) && nearly(tb[1][1], 0.0) && nearly(tb[1][2], 1.0));
    ASSERT_TRUE(nearly(tb[2][0], 0.0) && nearly(tb[2][1], -1.0) && nearly(tb[2][2], 0.0));

    /* HFDL: DICOM X=-uZ, DICOM Y=-uX, DICOM Z=+uY */
    osh_patient_position_base_rotation(OSH_PP_HFDL, tb);
    ASSERT_TRUE(nearly(tb[0][0], 0.0) && nearly(tb[0][1], 0.0) && nearly(tb[0][2], -1.0));
    ASSERT_TRUE(nearly(tb[1][0], -1.0) && nearly(tb[1][1], 0.0) && nearly(tb[1][2], 0.0));
    ASSERT_TRUE(nearly(tb[2][0], 0.0) && nearly(tb[2][1], 1.0) && nearly(tb[2][2], 0.0));

    /* HFDR: DICOM X=+uZ, DICOM Y=+uX, DICOM Z=+uY */
    osh_patient_position_base_rotation(OSH_PP_HFDR, tb);
    ASSERT_TRUE(nearly(tb[0][0], 0.0) && nearly(tb[0][1], 0.0) && nearly(tb[0][2], 1.0));
    ASSERT_TRUE(nearly(tb[1][0], 1.0) && nearly(tb[1][1], 0.0) && nearly(tb[1][2], 0.0));
    ASSERT_TRUE(nearly(tb[2][0], 0.0) && nearly(tb[2][1], 1.0) && nearly(tb[2][2], 0.0));

    /* FFDL: DICOM X=-uZ, DICOM Y=+uX, DICOM Z=-uY */
    osh_patient_position_base_rotation(OSH_PP_FFDL, tb);
    ASSERT_TRUE(nearly(tb[0][0], 0.0) && nearly(tb[0][1], 0.0) && nearly(tb[0][2], -1.0));
    ASSERT_TRUE(nearly(tb[1][0], 1.0) && nearly(tb[1][1], 0.0) && nearly(tb[1][2], 0.0));
    ASSERT_TRUE(nearly(tb[2][0], 0.0) && nearly(tb[2][1], -1.0) && nearly(tb[2][2], 0.0));

    /* FFDR: DICOM X=+uZ, DICOM Y=-uX, DICOM Z=-uY */
    osh_patient_position_base_rotation(OSH_PP_FFDR, tb);
    ASSERT_TRUE(nearly(tb[0][0], 0.0) && nearly(tb[0][1], 0.0) && nearly(tb[0][2], 1.0));
    ASSERT_TRUE(nearly(tb[1][0], -1.0) && nearly(tb[1][1], 0.0) && nearly(tb[1][2], 0.0));
    ASSERT_TRUE(nearly(tb[2][0], 0.0) && nearly(tb[2][1], -1.0) && nearly(tb[2][2], 0.0));
}

static void test_base_rotation_orthonormal(void) {
    enum osh_patient_position positions[] = {
        OSH_PP_HFS, OSH_PP_HFP, OSH_PP_FFS, OSH_PP_FFP, OSH_PP_HFDL, OSH_PP_HFDR, OSH_PP_FFDL, OSH_PP_FFDR};

    size_t i;
    double tb[3][3];

    for (i = 0; i < sizeof(positions) / sizeof(positions[0]); i++) {
        osh_patient_position_base_rotation(positions[i], tb);
        assert_proper_rotation(tb);
    }
}

static void test_base_rotation_default_fallback(void) {
    double tb_hfs[3][3];
    double tb_fallback[3][3];
    int i;
    int j;

    osh_patient_position_base_rotation(OSH_PP_HFS, tb_hfs);
    /* OSH_PP_UNKNOWN triggers the default branch, which should match HFS */
    osh_patient_position_base_rotation(OSH_PP_UNKNOWN, tb_fallback);
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            ASSERT_TRUE(nearly(tb_fallback[i][j], tb_hfs[i][j]));
}
