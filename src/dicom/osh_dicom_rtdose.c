#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_abort.h"
#include "common/osh_diag.h"
#include "dicom/osh_dicom_parse.h"
#include "openshieldhit/dicom.h"

struct _collect {
    struct osh_dicom_rtdose *rd;
    unsigned char const *buf; /* kept to compute pixel offset */
    struct osh_diag_sink const *diag;
    int ok;
};

static int _ds_value_capacity(unsigned char const *value, uint32_t length) {
    int count;
    uint32_t i;
    int in_token;

    if (!value || length == 0)
        return 0;

    count = 0;
    in_token = 0;
    for (i = 0; i < length; i++) {
        if (value[i] == '\\') {
            in_token = 0;
            continue;
        }
        if (!in_token && value[i] != ' ' && value[i] != '\0') {
            count++;
            in_token = 1;
        }
    }
    return count;
}

static int _rtdose_pixel_count(size_t *out, int n_frames, int rows, int cols) {
    size_t a;
    size_t b;
    size_t c;

    if (!out || n_frames <= 0 || rows <= 0 || cols <= 0)
        return 0;

    a = (size_t) n_frames;
    b = (size_t) rows;
    c = (size_t) cols;

    if (a > SIZE_MAX / b)
        return 0;
    a *= b;
    if (a > SIZE_MAX / c)
        return 0;
    *out = a * c;
    return 1;
}

static int
_tag_cb(uint16_t group, uint16_t element, char const vr[2], unsigned char const *value, uint32_t length, void *user) {
    struct _collect *c = (struct _collect *) user;
    struct osh_dicom_rtdose *rd = c->rd;
    (void) vr;

    if (group == 0x0020 && element == 0x0032) { /* Image Position Patient */
        double pos[3] = {0.0, 0.0, 0.0};
        osh_dicom_ds_array(value, length, pos, 3);
        rd->origin[0] = pos[0];
        rd->origin[1] = pos[1];
        rd->origin[2] = pos[2];
    } else if (group == 0x0020 && element == 0x0037) { /* Image Orientation Patient */
        double iop[6] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0};
        osh_dicom_ds_array(value, length, iop, 6);
        rd->row_cosine[0] = iop[0];
        rd->row_cosine[1] = iop[1];
        rd->row_cosine[2] = iop[2];
        rd->col_cosine[0] = iop[3];
        rd->col_cosine[1] = iop[4];
        rd->col_cosine[2] = iop[5];
    } else if (group == 0x0028 && element == 0x0008) { /* Number of Frames */
        char tmp[16] = {0};
        size_t n = length < 15 ? length : 15;
        memcpy(tmp, value, n);
        rd->n_frames = atoi(tmp);
    } else if (group == 0x0028 && element == 0x0010) { /* Rows */
        rd->rows = (int) osh_dicom_us(value, length);
    } else if (group == 0x0028 && element == 0x0011) { /* Columns */
        rd->cols = (int) osh_dicom_us(value, length);
    } else if (group == 0x0028 && element == 0x0030) { /* Pixel Spacing */
        osh_dicom_ds_array(value, length, rd->pixel_spacing, 2);
    } else if (group == 0x3004 && element == 0x000C) { /* Grid Frame Offset Vector */
        int capacity = _ds_value_capacity(value, length);
        int parsed;

        if (rd->n_frames > capacity)
            capacity = rd->n_frames;
        if (capacity < 1)
            capacity = 1;

        rd->frame_offsets = (double *) calloc((size_t) capacity, sizeof(double));
        if (!rd->frame_offsets)
            osh_abort_oomf("dicom rtdose: frame offsets");
        parsed = osh_dicom_ds_array(value, length, rd->frame_offsets, capacity);
        if (rd->n_frames > 0 && parsed != rd->n_frames) {
            OSH_DIAG_WARNF(
                c->diag, "dicom rtdose: parsed %d frame offsets but Number of Frames is %d", parsed, rd->n_frames);
        }
    } else if (group == 0x3004 && element == 0x000E) { /* Dose Grid Scaling */
        rd->dose_grid_scaling = osh_dicom_ds(value, length);
    } else if (group == 0x7FE0 && element == 0x0010) { /* Pixel Data */
        size_t offset;
        unsigned char *pixel_base;

        rd->n_pixels = 0;
        if (_rtdose_pixel_count(&rd->n_pixels, rd->n_frames, rd->rows, rd->cols)
            && length >= rd->n_pixels * sizeof(uint32_t)) {
            rd->_pixel_data_offset = (size_t) (value - c->buf);
            offset = rd->_pixel_data_offset;
            pixel_base = rd->_raw + offset;
            if (((uintptr_t) pixel_base % sizeof(uint32_t)) == 0u) {
                rd->pixels = (uint32_t *) pixel_base;
                c->ok = 1;
            }
        }
        return 0; /* always last — stop */
    }
    return 1;
}

enum osh_status osh_dicom_rtdose_read(char const *path, struct osh_dicom_rtdose *rd, struct osh_diag_sink const *diag) {
    unsigned char *buf;
    size_t buf_size;
    struct _collect c;

    memset(rd, 0, sizeof(*rd));

    buf = osh_dicom_load_file(path, &buf_size, diag);
    if (!buf)
        return OSH_EIO;

    rd->_raw = buf;
    rd->_raw_size = buf_size;

    c.rd = rd;
    c.buf = buf;
    c.diag = diag;
    c.ok = 0;

    osh_dicom_walk(buf, buf_size, _tag_cb, &c, diag);

    if (!c.ok) {
        OSH_DIAG_ERRORF(diag, "dicom rtdose: failed to read aligned pixel data from '%s'", path);
        osh_dicom_rtdose_free(rd);
        return OSH_EPARSE;
    }

    OSH_DIAG_INFOF(diag, "dicom rtdose: loaded %dx%dx%d from '%s'", rd->cols, rd->rows, rd->n_frames, path);
    return OSH_OK;
}

enum osh_status
osh_dicom_rtdose_write(char const *path, struct osh_dicom_rtdose const *rd, struct osh_diag_sink const *diag) {
    FILE *fp;

    if (!rd || !rd->_raw || rd->_raw_size == 0) {
        OSH_DIAG_ERRORF(diag, "dicom rtdose: nothing to write");
        return OSH_EINVAL;
    }

    fp = fopen(path, "wb");
    if (!fp) {
        OSH_DIAG_ERRORF(diag, "dicom rtdose: cannot open '%s' for writing", path);
        return OSH_EIO;
    }

    if (fwrite(rd->_raw, 1, rd->_raw_size, fp) != rd->_raw_size) {
        OSH_DIAG_ERRORF(diag, "dicom rtdose: write error on '%s'", path);
        fclose(fp);
        return OSH_EIO;
    }

    fclose(fp);
    OSH_DIAG_INFOF(diag, "dicom rtdose: wrote '%s'", path);
    return OSH_OK;
}

void osh_dicom_rtdose_free(struct osh_dicom_rtdose *rd) {
    if (!rd)
        return;
    free(rd->_raw); /* also invalidates rd->pixels */
    free(rd->frame_offsets);
    rd->_raw = NULL;
    rd->pixels = NULL;
    rd->frame_offsets = NULL;
}
