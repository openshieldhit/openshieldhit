#include <stdio.h>
#include <string.h>

#include "openshieldhit/diag.h"
#include "openshieldhit/dicom.h"

/*
 * Small developer utility for smoke-testing the DICOM reader.
 *
 * Keep this intentionally narrow: one local diagnostics sink plus two summary
 * commands (`ct` and `rtdose`). If the DICOM layer grows, document richer
 * workflows in README.md rather than turning this file into a framework.
 */
static void emit(void *user, int level, char const *file, int line, char const *function, char const *message) {
    (void) user;
    (void) file;
    (void) line;
    (void) function;
    fprintf(level >= OSH_DIAG_LEVEL_WARN ? stderr : stdout, "[%s] %s\n", osh_diag_level_name(level), message);
}

static struct osh_diag_sink diag = {emit, NULL, OSH_DIAG_LEVEL_INFO};

static void print_ct(char const *dir);
static void print_rtdose(char const *path);
static void usage(char const *prog);

int main(int argc, char *argv[]) {
    if (argc != 3) {
        usage(argv[0]);
        return 1;
    }

    if (strcmp(argv[1], "ct") == 0)
        print_ct(argv[2]);
    else if (strcmp(argv[1], "rtdose") == 0)
        print_rtdose(argv[2]);
    else {
        usage(argv[0]);
        return 1;
    }
    return 0;
}

static void usage(char const *prog) {
    fprintf(stderr, "usage: %s ct <directory>\n", prog);
    fprintf(stderr, "       %s rtdose <file>\n", prog);
}

static void print_ct(char const *dir) {
    struct osh_dicom_ct ct;
    int i;
    int16_t pmin, pmax;

    if (osh_dicom_ct_read(dir, &ct, &diag) != OSH_OK)
        return;

    printf("CT Series\n");
    printf("---------\n");
    printf("Directory     : %s\n", dir);
    printf("Dimensions    : %d cols x %d rows x %d slices\n", ct.cols, ct.rows, ct.n_slices);
    printf("Origin        : %9.3f  %9.3f  %9.3f mm\n", ct.origin[0], ct.origin[1], ct.origin[2]);
    printf("Row cosine    : %7.4f  %7.4f  %7.4f\n", ct.row_cosine[0], ct.row_cosine[1], ct.row_cosine[2]);
    printf("Col cosine    : %7.4f  %7.4f  %7.4f\n", ct.col_cosine[0], ct.col_cosine[1], ct.col_cosine[2]);
    printf("Pixel spacing : %.4f x %.4f mm (row x col)\n", ct.pixel_spacing[0], ct.pixel_spacing[1]);
    printf("Slice spacing : %.4f mm\n", ct.slice_spacing);
    printf("Rescale slope : %.6f\n", ct.rescale_slope);
    printf("Rescale intcp : %.6f\n", ct.rescale_intercept);

    pmin = pmax = ct.pixels[0];
    for (i = 1; i < ct.n_slices * ct.rows * ct.cols; i++) {
        if (ct.pixels[i] < pmin)
            pmin = ct.pixels[i];
        if (ct.pixels[i] > pmax)
            pmax = ct.pixels[i];
    }
    printf("Pixel range   : [%d, %d] (raw int16)\n", (int) pmin, (int) pmax);
    printf("HU range      : [%.0f, %.0f]\n",
           pmin * ct.rescale_slope + ct.rescale_intercept,
           pmax * ct.rescale_slope + ct.rescale_intercept);

    osh_dicom_ct_free(&ct);
}

static void print_rtdose(char const *path) {
    struct osh_dicom_rtdose rd;
    size_t i;
    uint32_t pmax;
    int fn;

    if (osh_dicom_rtdose_read(path, &rd, &diag) != OSH_OK)
        return;

    printf("RTDOSE\n");
    printf("------\n");
    printf("File          : %s\n", path);
    printf("Dimensions    : %d cols x %d rows x %d frames\n", rd.cols, rd.rows, rd.n_frames);
    printf("Origin        : %9.3f  %9.3f  %9.3f mm\n", rd.origin[0], rd.origin[1], rd.origin[2]);
    printf("Row cosine    : %7.4f  %7.4f  %7.4f\n", rd.row_cosine[0], rd.row_cosine[1], rd.row_cosine[2]);
    printf("Col cosine    : %7.4f  %7.4f  %7.4f\n", rd.col_cosine[0], rd.col_cosine[1], rd.col_cosine[2]);
    printf("Pixel spacing : %.4f x %.4f mm (row x col)\n", rd.pixel_spacing[0], rd.pixel_spacing[1]);
    printf("Grid scaling  : %.6e Gy/unit\n", rd.dose_grid_scaling);

    if (rd.frame_offsets) {
        printf("Frame offsets :");
        for (fn = 0; fn < rd.n_frames && fn < 8; fn++)
            printf(" %.2f", rd.frame_offsets[fn]);
        if (rd.n_frames > 8)
            printf(" ...");
        printf(" mm\n");
    }

    pmax = 0;
    for (i = 0; i < rd.n_pixels; i++)
        if (rd.pixels[i] > pmax)
            pmax = rd.pixels[i];
    printf("Max dose      : %.4f Gy\n", pmax * rd.dose_grid_scaling);

    osh_dicom_rtdose_free(&rd);
}
