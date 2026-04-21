#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_abort.h"
#include "common/osh_diag.h"
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

static int _slice_pixel_count(size_t *out, int rows, int cols) {
    size_t r;
    size_t c;

    if (!out || rows <= 0 || cols <= 0)
        return 0;

    r = (size_t) rows;
    c = (size_t) cols;
    if (r > SIZE_MAX / c)
        return 0;
    *out = r * c;
    return 1;
}

static int
_tag_cb(uint16_t group, uint16_t element, char const vr[2], unsigned char const *value, uint32_t length, void *user) {
    struct _collect *c = (struct _collect *) user;
    struct _slice *s = &c->s;
    (void) vr;

    if (group == 0x0008 && element == 0x0060) { /* Modality */
        char mod[17] = {0};
        size_t n = length < 16 ? length : 16;
        memcpy(mod, value, n);
        /* trim trailing spaces */
        while (n > 0 && mod[n - 1] == ' ')
            mod[--n] = '\0';
        s->is_ct = (strcmp(mod, "CT") == 0);
        if (!s->is_ct)
            return 0;                                  /* not CT — stop immediately */
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
        if (_slice_pixel_count(&n, s->rows, s->cols) && n <= (size_t) length / sizeof(int16_t)) {
            bytes = n * sizeof(int16_t);
            s->pixels = (int16_t *) malloc(bytes);
            if (!s->pixels)
                osh_abort_oomf("dicom ct: slice pixel buffer");
            memcpy(s->pixels, value, bytes);
        }
        return 0; /* pixel data is always last — stop walking */
    }
    return 1;
}

static int _z_cmp(void const *a, void const *b) {
    double za = ((struct _slice const *) a)->z;
    double zb = ((struct _slice const *) b)->z;
    return (za > zb) - (za < zb);
}

static int _is_dcm(char const *name) {
    size_t n = strlen(name);
    if (n < 4)
        return 0;
    char const *ext = name + n - 4;
    return (ext[0] == '.' && (ext[1] == 'd' || ext[1] == 'D') && (ext[2] == 'c' || ext[2] == 'C')
            && (ext[3] == 'm' || ext[3] == 'M'));
}

enum osh_status osh_dicom_ct_read(char const *dir, struct osh_dicom_ct *ct, struct osh_diag_sink const *diag) {
    DIR *d;
    struct dirent *ent;
    struct _slice *slices = NULL;
    int n = 0;
    int cap = 0;
    char path[4096];

    memset(ct, 0, sizeof(*ct));
    ct->rescale_slope = 1.0;

    d = opendir(dir);
    if (!d) {
        OSH_DIAG_ERRORF(diag, "dicom ct: cannot open directory '%s'", dir);
        return OSH_EIO;
    }

    while ((ent = readdir(d)) != NULL) {
        unsigned char *buf;
        size_t buf_size;
        struct _collect c;

        if (!_is_dcm(ent->d_name))
            continue;
        snprintf(path, sizeof(path), "%s/%s", dir, ent->d_name);

        buf = osh_dicom_load_file(path, &buf_size, diag);
        if (!buf)
            continue;

        memset(&c, 0, sizeof(c));
        c.s.rescale_slope = 1.0;
        osh_dicom_walk(buf, buf_size, _tag_cb, &c, diag);
        free(buf);

        if (!c.s.is_ct || !c.s.pixels || c.s.rows == 0 || c.s.cols == 0) {
            free(c.s.pixels);
            continue;
        }

        if (n == cap) {
            int new_cap = cap ? cap * 2 : 32;
            struct _slice *tmp = (struct _slice *) realloc(slices, (size_t) new_cap * sizeof(*slices));
            if (!tmp)
                osh_abort_oomf("dicom ct: slice array");
            slices = tmp;
            cap = new_cap;
        }
        slices[n++] = c.s;
    }
    closedir(d);

    if (n == 0) {
        OSH_DIAG_ERRORF(diag, "dicom ct: no valid CT slices found in '%s'", dir);
        free(slices);
        return OSH_EIO;
    }

    qsort(slices, (size_t) n, sizeof(*slices), _z_cmp);

    {
        int rows = slices[0].rows;
        int cols = slices[0].cols;
        size_t slice_px = (size_t) (rows * cols);
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
        if (!ct->pixels)
            osh_abort_oomf("dicom ct: pixel array");

        for (i = 0; i < n; i++) {
            memcpy(ct->pixels + (size_t) i * slice_px, slices[i].pixels, slice_px * sizeof(int16_t));
            free(slices[i].pixels);
        }
    }
    free(slices);

    OSH_DIAG_INFOF(diag, "dicom ct: loaded %d slices (%dx%d) from '%s'", n, ct->rows, ct->cols, dir);
    return OSH_OK;
}

void osh_dicom_ct_free(struct osh_dicom_ct *ct) {
    if (!ct)
        return;
    free(ct->pixels);
    ct->pixels = NULL;
}
