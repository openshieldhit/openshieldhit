#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_file.h"
#include "openshieldhit/dicom.h"
#include "openshieldhit/status.h"
#include "test_assert.h"

#define CT_DIR OSH_TEST_FIXTURES_DIR "/dicom/DCPT_headphantom"
#define RD_FILE                                                                                                        \
    OSH_TEST_FIXTURES_DIR "/dicom/DCPT_headphantom/"                                                                   \
                          "RD.1.2.246.352.71.7.37402163639.2864312.20240227185700.dcm"

static void test_ct_read_succeeds(void);
static void test_ct_geometry(void);
static void test_ct_rejects_non_uniform_slice_dimensions(void);
static void test_rtdose_read_succeeds(void);
static void test_rtdose_geometry(void);
static void test_rtdose_write_roundtrip(void);
static void write_u16_le(FILE *fp, uint16_t value);
static void write_u32_le(FILE *fp, uint32_t value);
static void write_tag_header(FILE *fp, uint16_t group, uint16_t element, char const vr[2], uint32_t length);
static void write_tag_string(FILE *fp, uint16_t group, uint16_t element, char const vr[2], char const *value);
static void write_tag_us(FILE *fp, uint16_t group, uint16_t element, uint16_t value);
static void write_tag_pixels(FILE *fp, int rows, int cols);
static void write_minimal_ct_slice(char const *path, int rows, int cols, char const *position);
static void join_path(char *out, size_t cap, char const *dir, char const *name);

int main(void) {
    test_ct_read_succeeds();
    test_ct_geometry();
    test_ct_rejects_non_uniform_slice_dimensions();
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
    double dot; /* Squared norm of the row direction cosine. */

    ASSERT_TRUE(osh_dicom_ct_read(CT_DIR, &ct, NULL) == OSH_OK);
    ASSERT_TRUE(ct.n_slices == 177);
    ASSERT_TRUE(ct.pixel_spacing[0] > 0.0);
    ASSERT_TRUE(ct.pixel_spacing[1] > 0.0);
    ASSERT_TRUE(ct.slice_spacing > 0.0);
    ASSERT_TRUE(ct.rescale_slope > 0.0);
    /* orientation cosines must be unit-ish (dot-product of row with itself ≈ 1) */
    dot =
        ct.row_cosine[0] * ct.row_cosine[0] + ct.row_cosine[1] * ct.row_cosine[1] + ct.row_cosine[2] * ct.row_cosine[2];
    ASSERT_TRUE(dot > 0.99 && dot < 1.01);
    osh_dicom_ct_free(&ct);
}

static void test_ct_rejects_non_uniform_slice_dimensions(void) {
    struct osh_dicom_ct ct;
    char const *dir = "osh_dicom_ct_mismatch.tmp";
    char slice_a[256]; /* First CT slice path, z=0, larger in-plane raster. */
    char slice_b[256]; /* Second CT slice path, z=1, smaller in-plane raster. */

    ASSERT_TRUE(osh_path_ensure_dir(dir) == OSH_OK);
    join_path(slice_a, sizeof(slice_a), dir, "slice_a.dcm");
    join_path(slice_b, sizeof(slice_b), dir, "slice_b.dcm");

    write_minimal_ct_slice(slice_a, 4, 4, "0\\0\\0 ");
    write_minimal_ct_slice(slice_b, 2, 2, "0\\0\\1 ");

    ASSERT_TRUE(osh_dicom_ct_read(dir, &ct, NULL) == OSH_EPARSE);
    osh_dicom_ct_free(&ct);

    (void) remove(slice_a);
    (void) remove(slice_b);
    (void) remove(dir);
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
    char const *tmp = "osh_test_rtdose_roundtrip.dcm";
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
    (void) remove(tmp);
}

static void write_u16_le(FILE *fp, uint16_t value) {
    ASSERT_TRUE(fputc((int) (value & 0xffu), fp) != EOF);
    ASSERT_TRUE(fputc((int) ((value >> 8) & 0xffu), fp) != EOF);
}

static void write_u32_le(FILE *fp, uint32_t value) {
    write_u16_le(fp, (uint16_t) (value & 0xffffu));
    write_u16_le(fp, (uint16_t) ((value >> 16) & 0xffffu));
}

static void write_tag_header(FILE *fp, uint16_t group, uint16_t element, char const vr[2], uint32_t length) {
    int is_long_vr; /* DICOM explicit-VR encodings with two reserved bytes and a 32-bit length. */

    is_long_vr = (vr[0] == 'O' && (vr[1] == 'B' || vr[1] == 'D' || vr[1] == 'F' || vr[1] == 'L' || vr[1] == 'W'))
                 || (vr[0] == 'S' && vr[1] == 'Q')
                 || (vr[0] == 'U' && (vr[1] == 'C' || vr[1] == 'R' || vr[1] == 'T' || vr[1] == 'N'));

    write_u16_le(fp, group);
    write_u16_le(fp, element);
    ASSERT_TRUE(fputc(vr[0], fp) != EOF);
    ASSERT_TRUE(fputc(vr[1], fp) != EOF);
    if (is_long_vr) {
        write_u16_le(fp, 0u);
        write_u32_le(fp, length);
    } else {
        ASSERT_TRUE(length <= UINT16_MAX);
        write_u16_le(fp, (uint16_t) length);
    }
}

static void write_tag_string(FILE *fp, uint16_t group, uint16_t element, char const vr[2], char const *value) {
    size_t len;        /* Raw byte count of the ASCII value before DICOM padding. */
    size_t padded_len; /* Even byte count written into the DICOM element header. */

    len = strlen(value);
    padded_len = len + (len % 2u);
    ASSERT_TRUE(padded_len <= UINT16_MAX);
    write_tag_header(fp, group, element, vr, (uint32_t) padded_len);
    ASSERT_TRUE(fwrite(value, 1u, len, fp) == len);
    if (padded_len != len) {
        ASSERT_TRUE(fputc(' ', fp) != EOF);
    }
}

static void write_tag_us(FILE *fp, uint16_t group, uint16_t element, uint16_t value) {
    write_tag_header(fp, group, element, "US", 2u);
    write_u16_le(fp, value);
}

static void write_tag_pixels(FILE *fp, int rows, int cols) {
    uint32_t n_pixels; /* Number of 16-bit pixels written for this synthetic slice. */
    uint32_t i;

    n_pixels = (uint32_t) rows * (uint32_t) cols;
    write_tag_header(fp, 0x7FE0u, 0x0010u, "OW", n_pixels * 2u);
    for (i = 0u; i < n_pixels; i++) {
        write_u16_le(fp, (uint16_t) i);
    }
}

static void write_minimal_ct_slice(char const *path, int rows, int cols, char const *position) {
    unsigned char preamble[128]; /* Standard DICOM preamble before the DICM marker. */
    FILE *fp;

    memset(preamble, 0, sizeof(preamble));
    fp = fopen(path, "wb");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fwrite(preamble, 1u, sizeof(preamble), fp) == sizeof(preamble));
    ASSERT_TRUE(fwrite("DICM", 1u, 4u, fp) == 4u);

    write_tag_string(fp, 0x0008u, 0x0060u, "CS", "CT");
    write_tag_string(fp, 0x0020u, 0x0032u, "DS", position);
    write_tag_string(fp, 0x0020u, 0x0037u, "DS", "1\\0\\0\\0\\1\\0 ");
    write_tag_us(fp, 0x0028u, 0x0010u, (uint16_t) rows);
    write_tag_us(fp, 0x0028u, 0x0011u, (uint16_t) cols);
    write_tag_string(fp, 0x0028u, 0x0030u, "DS", "1\\1 ");
    write_tag_string(fp, 0x0028u, 0x1052u, "DS", "0 ");
    write_tag_string(fp, 0x0028u, 0x1053u, "DS", "1 ");
    write_tag_pixels(fp, rows, cols);

    ASSERT_TRUE(fclose(fp) == 0);
}

static void join_path(char *out, size_t cap, char const *dir, char const *name) {
    int n; /* snprintf result used to reject truncated test paths. */

    n = snprintf(out, cap, "%s/%s", dir, name);
    ASSERT_TRUE(n >= 0 && n < (int) cap);
}
