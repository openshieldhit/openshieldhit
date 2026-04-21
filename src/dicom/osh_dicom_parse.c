#include "dicom/osh_dicom_parse.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_abort.h"
#include "common/osh_diag.h"

static uint16_t _u16(unsigned char const *p) {
    return (uint16_t) ((unsigned) p[0] | ((unsigned) p[1] << 8));
}

static uint32_t _u32(unsigned char const *p) {
    return (uint32_t) ((unsigned) p[0] | ((unsigned) p[1] << 8) | ((unsigned) p[2] << 16) | ((unsigned) p[3] << 24));
}

static int _is_long_vr(unsigned char const vr[2]) {
    /* VRs that encode length as 4 bytes (with 2 reserved bytes): */
    /* OB OD OF OL OV OW SQ UC UN UR UT */
    switch (vr[0]) {
    case 'O':
        return vr[1] == 'B' || vr[1] == 'D' || vr[1] == 'F' || vr[1] == 'L' || vr[1] == 'V' || vr[1] == 'W';
    case 'S':
        return vr[1] == 'Q';
    case 'U':
        return vr[1] == 'C' || vr[1] == 'N' || vr[1] == 'R' || vr[1] == 'T';
    default:
        return 0;
    }
}

unsigned char *osh_dicom_load_file(char const *path, size_t *out_size, struct osh_diag_sink const *diag) {
    FILE *fp;
    long sz;
    unsigned char *buf;

    fp = fopen(path, "rb");
    if (!fp) {
        OSH_DIAG_ERRORF(diag, "dicom: cannot open '%s'", path);
        return NULL;
    }
    if (fseek(fp, 0, SEEK_END) != 0 || (sz = ftell(fp)) < 0) {
        OSH_DIAG_ERRORF(diag, "dicom: cannot determine size of '%s'", path);
        fclose(fp);
        return NULL;
    }
    rewind(fp);

    buf = (unsigned char *) malloc((size_t) sz);
    if (!buf)
        osh_abort_oomf("dicom: loading '%s'", path);

    if (fread(buf, 1, (size_t) sz, fp) != (size_t) sz) {
        OSH_DIAG_ERRORF(diag, "dicom: read error on '%s'", path);
        free(buf);
        fclose(fp);
        return NULL;
    }
    fclose(fp);
    *out_size = (size_t) sz;
    return buf;
}

int osh_dicom_walk(
    unsigned char const *buf, size_t size, osh_dicom_tag_fn fn, void *user, struct osh_diag_sink const *diag) {
    size_t pos;

    if (size < 132 || memcmp(buf + 128, "DICM", 4) != 0) {
        OSH_DIAG_ERRORF(diag, "dicom: missing DICM preamble");
        return 0;
    }
    pos = 132;

    while (pos + 6 <= size) {
        uint16_t group = _u16(buf + pos);
        uint16_t element = _u16(buf + pos + 2);
        unsigned char const *vr = buf + pos + 4;
        uint32_t length;
        size_t value_start;

        pos += 6; /* group(2) + element(2) + vr(2) */

        if (_is_long_vr(vr)) {
            if (pos + 6 > size) {
                OSH_DIAG_ERRORF(
                    diag, "dicom: truncated long-VR header at (%04X,%04X)", (unsigned) group, (unsigned) element);
                return 0;
            }
            pos += 2; /* reserved bytes */
            length = _u32(buf + pos);
            pos += 4;
        } else {
            if (pos + 2 > size) {
                OSH_DIAG_ERRORF(
                    diag, "dicom: truncated short-VR header at (%04X,%04X)", (unsigned) group, (unsigned) element);
                return 0;
            }
            length = _u16(buf + pos);
            pos += 2;
        }

        /* Undefined-length items (sequences) — not needed for flat tags */
        if (length == 0xFFFFFFFFu) {
            OSH_DIAG_WARNF(diag, "dicom: undefined-length item at (%04X,%04X), stopping walk", group, element);
            break;
        }

        if (pos + length > size) {
            OSH_DIAG_ERRORF(diag,
                            "dicom: tag (%04X,%04X) overruns buffer (%u bytes at offset %zu)",
                            (unsigned) group,
                            (unsigned) element,
                            (unsigned) length,
                            pos);
            return 0;
        }
        value_start = pos;
        pos += length;

        if (!fn(group, element, (char const *) vr, buf + value_start, length, user))
            return 1; /* caller requested stop */
    }
    return 1;
}

int osh_dicom_ds_array(unsigned char const *value, uint32_t length, double *out, int max_count) {
    char tmp[64];
    int count = 0;
    uint32_t i = 0;

    while (i < length && count < max_count) {
        size_t j = 0;
        /* skip leading whitespace / null padding */
        while (i < length && (value[i] == ' ' || value[i] == '\0'))
            i++;
        if (i >= length)
            break;
        /* copy token up to backslash separator */
        while (i < length && value[i] != '\\' && j < sizeof(tmp) - 1)
            tmp[j++] = (char) value[i++];
        tmp[j] = '\0';
        if (j > 0)
            out[count++] = atof(tmp);
        if (i < length && value[i] == '\\')
            i++;
    }
    return count;
}

double osh_dicom_ds(unsigned char const *value, uint32_t length) {
    double v = 0.0;
    osh_dicom_ds_array(value, length, &v, 1);
    return v;
}

uint16_t osh_dicom_us(unsigned char const *value, uint32_t length) {
    if (length < 2)
        return 0;
    return _u16(value);
}
