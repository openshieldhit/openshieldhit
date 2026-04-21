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
    int ok;
};

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
        int count = length / 12;                       /* rough upper bound for DS array */
        if (count < 1)
            count = rd->n_frames > 0 ? rd->n_frames : 256;
        rd->frame_offsets = (double *) malloc((size_t) count * sizeof(double));
        if (!rd->frame_offsets)
            osh_abort_oomf("dicom rtdose: frame offsets");
        osh_dicom_ds_array(value, length, rd->frame_offsets, count);
    } else if (group == 0x3004 && element == 0x000E) { /* Dose Grid Scaling */
        rd->dose_grid_scaling = osh_dicom_ds(value, length);
    } else if (group == 0x7FE0 && element == 0x0010) { /* Pixel Data */
        rd->n_pixels = (size_t) (rd->n_frames * rd->rows * rd->cols);
        if (rd->n_pixels > 0 && length >= rd->n_pixels * sizeof(uint32_t)) {
            rd->_pixel_data_offset = (size_t) (value - c->buf);
            rd->pixels = (uint32_t *) (rd->_raw + rd->_pixel_data_offset);
            c->ok = 1;
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
    c.ok = 0;

    osh_dicom_walk(buf, buf_size, _tag_cb, &c, diag);

    if (!c.ok) {
        OSH_DIAG_ERRORF(diag, "dicom rtdose: failed to read pixel data from '%s'", path);
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
