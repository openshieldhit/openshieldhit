#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_abort.h"
#include "common/osh_diag.h"
#include "common/osh_file.h"
#include "dicom/osh_dicom_parse.h"
#include "openshieldhit/dicom.h"

/* Per-slice data collected from one DICOM file */
struct _slice {
    double z;
    double origin[3];
    double row_cosine[3];
    double col_cosine[3];
    double pixel_spacing[2];
    double rescale_intercept;
    double rescale_slope;
    int rows;
    int cols;
    int is_ct; /* modality == "CT" */
    int16_t *pixels;
};

struct _collect {
    struct _slice s;
};

struct _scan_ct_dir {
    struct _slice **slices;
    int *n;
    int *cap;
    enum osh_status rc;
    struct osh_diag_sink const *diag;
};

/**
 * @brief Free the temporary per-file CT slices collected while scanning a directory.
 *
 * @details
 * Each accepted file owns one pixel buffer until the series is validated and
 * copied into the public contiguous CT image.  Error exits before that copy
 * must release every still-owned slice buffer.
 */
static void _free_slices(struct _slice *slices, int n) {
    int i;

    if (!slices) {
        return;
    }
    for (i = 0; i < n; i++) {
        free(slices[i].pixels);
    }
    free(slices);
}

/**
 * @brief Compute the number of pixels in one CT slice.
 *
 * @details
 * Rows and columns come directly from DICOM unsigned-short tags, but are stored
 * as signed ints in the temporary slice state so missing or invalid values can
 * be represented as non-positive.  This helper rejects those values and guards
 * the rows*cols multiplication before callers use the result for allocation or
 * copy bounds.
 *
 * @param[out] out   Pixel count for one slice on success.
 * @param[in]  rows  DICOM Rows value.
 * @param[in]  cols  DICOM Columns value.
 *
 * @returns OSH_OK on success, OSH_EINVAL for invalid dimensions or overflow.
 */
static enum osh_status _slice_pixel_count(size_t *out, int rows, int cols) {
    size_t r;
    size_t c;

    if (!out || rows <= 0 || cols <= 0) {
        return OSH_EINVAL;
    }
    r = (size_t) rows;
    c = (size_t) cols;
    if (r > SIZE_MAX / c) {
        return OSH_EINVAL;
    }
    *out = r * c;
    return OSH_OK;
}

/**
 * @brief Collect the CT tags needed from one DICOM file.
 *
 * @details
 * This is an osh_dicom_walk() callback, not a status-returning operation.
 * Returning 1 asks the walker to continue; returning 0 stops walking the file.
 * Non-CT files stop immediately after Modality, and CT files stop after Pixel
 * Data because the reader only needs the flat tags before and including pixels.
 *
 * @returns 1 to continue DICOM tag walking, 0 to stop walking this file.
 */
static int
_tag_cb(uint16_t group, uint16_t element, char const vr[2], unsigned char const *value, uint32_t length, void *user) {
    struct _collect *c = (struct _collect *) user;
    struct _slice *s = &c->s;
    (void) vr;

    if (group == 0x0008 && element == 0x0060) { /* Modality */
        char mod[17] = {0};
        size_t n = length < 16 ? length : 16;
        memcpy(mod, value, n);
        while (n > 0 && mod[n - 1] == ' ') {
            mod[--n] = '\0';
        }
        s->is_ct = (strcmp(mod, "CT") == 0);
        if (!s->is_ct) {
            return 0; /* not CT — stop immediately */
        }
    } else if (group == 0x0020 && element == 0x0032) { /* Image Position Patient */
        double pos[3] = {0.0, 0.0, 0.0};
        osh_dicom_ds_array(value, length, pos, 3);
        s->origin[0] = pos[0];
        s->origin[1] = pos[1];
        s->origin[2] = pos[2];
        s->z = pos[2];
    } else if (group == 0x0020 && element == 0x0037) { /* Image Orientation Patient */
        double iop[6] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0};
        osh_dicom_ds_array(value, length, iop, 6);
        s->row_cosine[0] = iop[0];
        s->row_cosine[1] = iop[1];
        s->row_cosine[2] = iop[2];
        s->col_cosine[0] = iop[3];
        s->col_cosine[1] = iop[4];
        s->col_cosine[2] = iop[5];
    } else if (group == 0x0028 && element == 0x0010) { /* Rows */
        s->rows = (int) osh_dicom_us(value, length);
    } else if (group == 0x0028 && element == 0x0011) { /* Columns */
        s->cols = (int) osh_dicom_us(value, length);
    } else if (group == 0x0028 && element == 0x0030) { /* Pixel Spacing */
        osh_dicom_ds_array(value, length, s->pixel_spacing, 2);
    } else if (group == 0x0028 && element == 0x1052) { /* Rescale Intercept */
        s->rescale_intercept = osh_dicom_ds(value, length);
    } else if (group == 0x0028 && element == 0x1053) { /* Rescale Slope */
        s->rescale_slope = osh_dicom_ds(value, length);
    } else if (group == 0x7FE0 && element == 0x0010) { /* Pixel Data */
        size_t n;
        size_t bytes;
        if (_slice_pixel_count(&n, s->rows, s->cols) == OSH_OK && n <= (size_t) length / sizeof(int16_t)) {
            bytes = n * sizeof(int16_t);
            s->pixels = (int16_t *) malloc(bytes);
            if (!s->pixels) {
                osh_abort_oomf("dicom ct: slice pixel buffer");
            }
            memcpy(s->pixels, value, bytes);
        }
        return 0; /* pixel data is always last — stop walking */
    }
    return 1;
}

/**
 * @brief qsort() comparator for ascending CT slice z position.
 *
 * @returns Negative, zero, or positive according to the qsort comparator
 *          contract.
 */
static int _z_cmp(void const *a, void const *b) {
    double za = ((struct _slice const *) a)->z;
    double zb = ((struct _slice const *) b)->z;
    return (za > zb) - (za < zb);
}

/**
 * @brief Return whether a path has a .dcm extension.
 *
 * @returns 1 for .dcm/.DCM-style names, 0 otherwise.
 */
static int _is_dcm(char const *name) {
    size_t n = strlen(name);
    if (n < 4) {
        return 0;
    }
    char const *ext = name + n - 4;
    return (ext[0] == '.' && (ext[1] == 'd' || ext[1] == 'D') && (ext[2] == 'c' || ext[2] == 'C')
            && (ext[3] == 'm' || ext[3] == 'M'));
}

/**
 * @brief Read one candidate file and append it if it is a valid CT slice.
 *
 * @details
 * Non-DICOM files, non-CT DICOM files, and CT files lacking usable Rows,
 * Columns, or Pixel Data are ignored.  Accepted slices transfer ownership of
 * their pixel buffer into @p slices; ignored slices free any partial pixel
 * buffer before returning.
 *
 * @returns OSH_OK. Allocation failures abort via osh_abort_oomf().
 */
static enum osh_status
_append_ct_slice(struct _slice **slices, int *n, int *cap, char const *path, struct osh_diag_sink const *diag) {
    unsigned char *buf;
    size_t buf_size;
    struct _collect c;

    buf = osh_dicom_load_file(path, &buf_size, diag);
    if (!buf) {
        return OSH_OK;
    }

    memset(&c, 0, sizeof(c));
    c.s.rescale_slope = 1.0;
    osh_dicom_walk(buf, buf_size, _tag_cb, &c, diag);
    free(buf);

    if (!c.s.is_ct || !c.s.pixels || c.s.rows == 0 || c.s.cols == 0) {
        free(c.s.pixels);
        return OSH_OK;
    }

    if (*n == *cap) {
        int new_cap = *cap ? *cap * 2 : 32;
        struct _slice *tmp = (struct _slice *) realloc(*slices, (size_t) new_cap * sizeof(**slices));
        if (!tmp) {
            osh_abort_oomf("dicom ct: slice array");
        }
        *slices = tmp;
        *cap = new_cap;
    }
    (*slices)[(*n)++] = c.s;
    return OSH_OK;
}

/**
 * @brief Directory-walker callback that scans .dcm files for CT slices.
 *
 * @details
 * This is an osh_dir_foreach_file() callback.  Returning 1 continues directory
 * iteration; returning 0 stops it.  Operation failures from _append_ct_slice()
 * are stored in the scan context so osh_dicom_ct_read() can return the status
 * after osh_dir_foreach_file() unwinds.
 *
 * @returns 1 to continue directory iteration, 0 to stop after an operation
 *          failure.
 */
static int _scan_ct_file(char const *path, void *user) {
    struct _scan_ct_dir *scan = (struct _scan_ct_dir *) user;
    enum osh_status rc;

    if (!_is_dcm(path)) {
        return 1;
    }
    rc = _append_ct_slice(scan->slices, scan->n, scan->cap, path, scan->diag);
    if (rc != OSH_OK) {
        scan->rc = rc;
        return 0;
    }
    return 1;
}

/**
 * @brief Reject CT series whose slices do not all share one in-plane raster.
 *
 * @details
 * osh_dicom_ct_read() stores the series as one flat rows*cols*n_slices array.
 * Copying a mixed-size series through that layout would either read past a
 * smaller slice buffer or silently truncate a larger one, so this invariant is
 * checked immediately after z-sorting and before the public CT image is built.
 */
static enum osh_status
_validate_uniform_dimensions(struct _slice const *slices, int n, char const *dir, struct osh_diag_sink const *diag) {
    int rows; /* Reference row count from the first z-sorted slice. */
    int cols; /* Reference column count from the first z-sorted slice. */
    int i;

    if (!slices || n <= 0) {
        return OSH_EINVAL;
    }

    rows = slices[0].rows;
    cols = slices[0].cols;
    for (i = 1; i < n; i++) {
        if (slices[i].rows == rows && slices[i].cols == cols) {
            continue;
        }
        OSH_DIAG_ERRORF(diag,
                        "dicom ct: inconsistent slice dimensions in '%s' "
                        "(slice 0 is %dx%d, slice %d is %dx%d)",
                        dir,
                        rows,
                        cols,
                        i,
                        slices[i].rows,
                        slices[i].cols);
        return OSH_EPARSE;
    }

    return OSH_OK;
}

enum osh_status osh_dicom_ct_read(char const *dir, struct osh_dicom_ct *ct, struct osh_diag_sink const *diag) {
    struct _slice *slices = NULL;
    int n = 0;
    int cap = 0;
    struct _scan_ct_dir scan;
    enum osh_status rc;

    memset(ct, 0, sizeof(*ct));
    ct->rescale_slope = 1.0;

    scan.slices = &slices;
    scan.n = &n;
    scan.cap = &cap;
    scan.rc = OSH_OK;
    scan.diag = diag;

    rc = osh_dir_foreach_file(dir, _scan_ct_file, &scan);
    if (rc != OSH_OK) {
        OSH_DIAG_ERRORF(diag, "dicom ct: cannot open directory '%s'", dir);
        _free_slices(slices, n);
        return OSH_EIO;
    }
    if (scan.rc != OSH_OK) {
        _free_slices(slices, n);
        return scan.rc;
    }

    if (n == 0) {
        OSH_DIAG_ERRORF(diag, "dicom ct: no valid CT slices found in '%s'", dir);
        _free_slices(slices, n);
        return OSH_EIO;
    }

    qsort(slices, (size_t) n, sizeof(*slices), _z_cmp);
    rc = _validate_uniform_dimensions(slices, n, dir, diag);
    if (rc != OSH_OK) {
        _free_slices(slices, n);
        return rc;
    }

    {
        int rows = slices[0].rows;
        int cols = slices[0].cols;
        size_t slice_px = (size_t) rows * (size_t) cols;
        int16_t *dst;
        int16_t const *src;
        double slope;
        double intercept;
        size_t j;
        int i;

        ct->rows = rows;
        ct->cols = cols;
        ct->n_slices = n;
        ct->origin[0] = slices[0].origin[0];
        ct->origin[1] = slices[0].origin[1];
        ct->origin[2] = slices[0].origin[2];
        ct->row_cosine[0] = slices[0].row_cosine[0];
        ct->row_cosine[1] = slices[0].row_cosine[1];
        ct->row_cosine[2] = slices[0].row_cosine[2];
        ct->col_cosine[0] = slices[0].col_cosine[0];
        ct->col_cosine[1] = slices[0].col_cosine[1];
        ct->col_cosine[2] = slices[0].col_cosine[2];
        ct->pixel_spacing[0] = slices[0].pixel_spacing[0];
        ct->pixel_spacing[1] = slices[0].pixel_spacing[1];
        ct->rescale_intercept = slices[0].rescale_intercept;
        ct->rescale_slope = slices[0].rescale_slope;
        ct->slice_spacing = (n > 1) ? (slices[n - 1].z - slices[0].z) / (double) (n - 1) : 0.0;

        ct->pixels = (int16_t *) malloc(slice_px * (size_t) n * sizeof(int16_t));
        if (!ct->pixels) {
            osh_abort_oomf("dicom ct: pixel array");
        }

        for (i = 0; i < n; i++) {
            dst = ct->pixels + ((size_t) i * slice_px);
            src = slices[i].pixels;
            slope = slices[i].rescale_slope;
            intercept = slices[i].rescale_intercept;
            for (j = 0; j < slice_px; j++) {
                double hu = (double) src[j] * slope + intercept;
                if (hu < -32768.0)
                    hu = -32768.0;
                if (hu > 32767.0)
                    hu = 32767.0;
                dst[j] = (int16_t) hu;
            }
            free(slices[i].pixels);
        }
    }
    free(slices);

    OSH_DIAG_INFOF(diag, "dicom ct: loaded %d slices (%dx%d) from '%s'", n, ct->rows, ct->cols, dir);
    return OSH_OK;
}

void osh_dicom_ct_free(struct osh_dicom_ct *ct) {
    if (!ct) {
        return;
    }
    free(ct->pixels);
    ct->pixels = NULL;
}
