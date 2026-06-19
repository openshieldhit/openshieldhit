#include "scoring/save/osh_scoring_save_rtdose.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_diag.h"
#include "openshieldhit/dicom.h"
#include "scoring/runtime/osh_scoring_geometry_runtime.h"
#include "scoring/runtime/osh_scoring_output_runtime.h"
#include "scoring/runtime/osh_scoring_runtime.h"

static char *rtdose_output_path(char const *filename);

enum osh_status osh_scoring_save_rtdose_output(struct osh_scoring_workspace const *ws,
                                               struct osh_scoring_runtime const *rt,
                                               unsigned long long nstat,
                                               size_t output_idx) {
    struct osh_scoring_output_runtime const *out;
    struct osh_scoring_geometry_runtime const *geo;
    struct osh_scoring_page_runtime const *page;
    struct osh_dicom_rtdose rd;
    char *out_path = NULL;
    size_t expected_bins;
    size_t i;
    double inv_nstat;
    double inv_scaling;
    enum osh_status rc;

    (void) ws;

    if (!rt || output_idx >= rt->noutputs) {
        return OSH_EINVAL;
    }
    if (nstat == 0ull) {
        return OSH_EINVAL;
    }

    out = &rt->outputs[output_idx];
    geo = &rt->geometries[out->geometry_idx];

    if (!geo->rtdose_template_path) {
        return OSH_EINVAL;
    }
    if (out->npages != 1u) {
        return OSH_ENOTSUP;
    }

    page = &rt->pages[out->page_indices[0]];

    rc = osh_dicom_rtdose_read(geo->rtdose_template_path, &rd, NULL);
    if (rc != OSH_OK) {
        return rc;
    }

    expected_bins = (size_t) rd.n_frames * (size_t) rd.rows * (size_t) rd.cols;
    if (expected_bins != geo->nbins || expected_bins != page->len) {
        osh_dicom_rtdose_free(&rd);
        return OSH_ESTATE;
    }

    inv_nstat = 1.0 / (double) nstat;
    inv_scaling = (rd.dose_grid_scaling > 0.0) ? (1.0 / rd.dose_grid_scaling) : 1.0;

    for (i = 0; i < page->len; ++i) {
        double val = page->acc.data[i] * inv_nstat * inv_scaling;
        if (val < 0.0) {
            val = 0.0;
        }
        if (val > 4294967295.0) {
            val = 4294967295.0;
        }
        rd.pixels[i] = (uint32_t) val;
    }

    out_path = rtdose_output_path(out->filename);
    if (!out_path) {
        osh_dicom_rtdose_free(&rd);
        return OSH_ENOMEM;
    }

    rc = osh_dicom_rtdose_write(out_path, &rd, NULL);
    free(out_path);
    osh_dicom_rtdose_free(&rd);
    return rc;
}

/* Append ".dcm" to filename if not already present. */
static char *rtdose_output_path(char const *filename) {
    size_t len;
    char *result;

    if (!filename) {
        return NULL;
    }
    len = strlen(filename);
    if (len >= 4u && strcmp(filename + len - 4u, ".dcm") == 0) {
        result = (char *) malloc(len + 1u);
        if (result) {
            memcpy(result, filename, len + 1u);
        }
        return result;
    }
    result = (char *) malloc(len + 5u);
    if (result) {
        memcpy(result, filename, len);
        memcpy(result + len, ".dcm", 5u);
    }
    return result;
}
