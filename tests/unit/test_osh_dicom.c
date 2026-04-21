#include <stddef.h>

#include "openshieldhit/dicom.h"
#include "openshieldhit/status.h"
#include "test_assert.h"

#define CT_DIR OSH_TEST_FIXTURES_DIR "/dicom/DCPT_headphantom"
#define RD_FILE                                                                                                        \
    OSH_TEST_FIXTURES_DIR "/dicom/DCPT_headphantom/"                                                                   \
                          "RD.1.2.246.352.71.7.37402163639.2864312.20240227185700.dcm"

static void test_ct_read_succeeds(void);
static void test_ct_geometry(void);
static void test_rtdose_read_succeeds(void);
static void test_rtdose_geometry(void);
static void test_rtdose_write_roundtrip(void);

int main(void) {
    test_ct_read_succeeds();
    test_ct_geometry();
    test_rtdose_read_succeeds();
    test_rtdose_geometry();
    test_rtdose_write_roundtrip();
    return 0;
}

static void test_ct_read_succeeds(void) {
    struct osh_dicom_ct ct;
    ASSERT_TRUE(osh_dicom_ct_read(CT_DIR, &ct, NULL) == OSH_OK);
    ASSERT_TRUE(ct.n_slices > 0);
    ASSERT_TRUE(ct.rows > 0);
    ASSERT_TRUE(ct.cols > 0);
    ASSERT_TRUE(ct.pixels != NULL);
    osh_dicom_ct_free(&ct);
}

static void test_ct_geometry(void) {
    struct osh_dicom_ct ct;
    ASSERT_TRUE(osh_dicom_ct_read(CT_DIR, &ct, NULL) == OSH_OK);
    ASSERT_TRUE(ct.n_slices == 177);
    ASSERT_TRUE(ct.pixel_spacing[0] > 0.0);
    ASSERT_TRUE(ct.pixel_spacing[1] > 0.0);
    ASSERT_TRUE(ct.slice_spacing > 0.0);
    ASSERT_TRUE(ct.rescale_slope > 0.0);
    /* orientation cosines must be unit-ish (dot-product of row with itself ≈ 1) */
    double dot =
        ct.row_cosine[0] * ct.row_cosine[0] + ct.row_cosine[1] * ct.row_cosine[1] + ct.row_cosine[2] * ct.row_cosine[2];
    ASSERT_TRUE(dot > 0.99 && dot < 1.01);
    osh_dicom_ct_free(&ct);
}

static void test_rtdose_read_succeeds(void) {
    struct osh_dicom_rtdose rd;
    ASSERT_TRUE(osh_dicom_rtdose_read(RD_FILE, &rd, NULL) == OSH_OK);
    ASSERT_TRUE(rd.n_frames > 0);
    ASSERT_TRUE(rd.rows > 0);
    ASSERT_TRUE(rd.cols > 0);
    ASSERT_TRUE(rd.pixels != NULL);
    ASSERT_TRUE(rd.n_pixels > 0);
    osh_dicom_rtdose_free(&rd);
}

static void test_rtdose_geometry(void) {
    struct osh_dicom_rtdose rd;
    ASSERT_TRUE(osh_dicom_rtdose_read(RD_FILE, &rd, NULL) == OSH_OK);
    ASSERT_TRUE(rd.dose_grid_scaling > 0.0);
    ASSERT_TRUE(rd.pixel_spacing[0] > 0.0);
    ASSERT_TRUE(rd.pixel_spacing[1] > 0.0);
    ASSERT_TRUE(rd.frame_offsets != NULL);
    osh_dicom_rtdose_free(&rd);
}

static void test_rtdose_write_roundtrip(void) {
    struct osh_dicom_rtdose rd;
    struct osh_dicom_rtdose rd2;
    char const *tmp = "/tmp/osh_test_rtdose_roundtrip.dcm";
    size_t i;

    ASSERT_TRUE(osh_dicom_rtdose_read(RD_FILE, &rd, NULL) == OSH_OK);
    ASSERT_TRUE(osh_dicom_rtdose_write(tmp, &rd, NULL) == OSH_OK);

    /* re-read and verify pixel values are identical */
    ASSERT_TRUE(osh_dicom_rtdose_read(tmp, &rd2, NULL) == OSH_OK);
    ASSERT_TRUE(rd2.n_pixels == rd.n_pixels);
    for (i = 0; i < rd.n_pixels; i++)
        ASSERT_TRUE(rd2.pixels[i] == rd.pixels[i]);

    osh_dicom_rtdose_free(&rd);
    osh_dicom_rtdose_free(&rd2);
}
