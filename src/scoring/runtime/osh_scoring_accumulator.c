#include "scoring/runtime/osh_scoring_accumulator.h"

#include <stdlib.h>
#include <string.h>

enum osh_status osh_scoring_accumulator_alloc(struct osh_scoring_accumulator *acc, size_t len, int want_data2) {
    size_t const n = len ? len : 1u; /* never allocate a zero-length array */

    if (!acc) {
        return OSH_EINVAL;
    }

    acc->data = NULL;
    acc->data2 = NULL;
    acc->data_var = NULL;
    acc->data2_var = NULL;
    acc->len = len;

    acc->data = (double *) calloc(n, sizeof(*acc->data));
    if (!acc->data) {
        acc->len = 0u;
        return OSH_ENOMEM;
    }
    if (want_data2) {
        acc->data2 = (double *) calloc(n, sizeof(*acc->data2));
        if (!acc->data2) {
            free(acc->data);
            acc->data = NULL;
            acc->len = 0u;
            return OSH_ENOMEM;
        }
    }
    return OSH_OK;
}

void osh_scoring_accumulator_zero(struct osh_scoring_accumulator *acc) {
    if (!acc || acc->len == 0u) {
        return;
    }
    if (acc->data) {
        memset(acc->data, 0, acc->len * sizeof(*acc->data));
    }
    if (acc->data2) {
        memset(acc->data2, 0, acc->len * sizeof(*acc->data2));
    }
    if (acc->data_var) {
        memset(acc->data_var, 0, acc->len * sizeof(*acc->data_var));
    }
    if (acc->data2_var) {
        memset(acc->data2_var, 0, acc->len * sizeof(*acc->data2_var));
    }
}

void osh_scoring_accumulator_rescale(struct osh_scoring_accumulator *acc, double factor) {
    size_t i;
    if (!acc || !acc->data) {
        return;
    }
    for (i = 0; i < acc->len; ++i) {
        acc->data[i] *= factor;
    }
}

enum osh_status osh_scoring_accumulator_finalize_average(struct osh_scoring_accumulator *acc, double eps) {
    size_t i;
    if (!acc || !acc->data || !acc->data2) {
        return OSH_EINVAL;
    }
    for (i = 0; i < acc->len; ++i) {
        acc->data[i] = (acc->data2[i] > eps) ? (acc->data[i] / acc->data2[i]) : 0.0;
    }
    return OSH_OK;
}

/* Add src_arr into dst_arr element-wise when both are non-NULL. */
static void merge_array(double *dst_arr, double const *src_arr, size_t len) {
    size_t i;
    if (!dst_arr || !src_arr) {
        return;
    }
    for (i = 0; i < len; ++i) {
        dst_arr[i] += src_arr[i];
    }
}

enum osh_status osh_scoring_accumulator_merge(struct osh_scoring_accumulator *dst,
                                              struct osh_scoring_accumulator const *src) {
    if (!dst || !src) {
        return OSH_EINVAL;
    }
    if (dst->len != src->len) {
        return OSH_EINVAL;
    }
    /* Reject mismatched optional-array presence: an array allocated on one side
     * but NULL on the other would silently drop tallies.  Accumulators cloned
     * from the same page descriptor always agree, so disagreement is a bug. */
    if (((dst->data == NULL) != (src->data == NULL)) || ((dst->data2 == NULL) != (src->data2 == NULL))
        || ((dst->data_var == NULL) != (src->data_var == NULL))
        || ((dst->data2_var == NULL) != (src->data2_var == NULL))) {
        return OSH_EINVAL;
    }
    merge_array(dst->data, src->data, dst->len);
    merge_array(dst->data2, src->data2, dst->len);
    merge_array(dst->data_var, src->data_var, dst->len);
    merge_array(dst->data2_var, src->data2_var, dst->len);
    return OSH_OK;
}

void osh_scoring_accumulator_free(struct osh_scoring_accumulator *acc) {
    if (!acc) {
        return;
    }
    free(acc->data);
    free(acc->data2);
    free(acc->data_var);
    free(acc->data2_var);
    acc->data = NULL;
    acc->data2 = NULL;
    acc->data_var = NULL;
    acc->data2_var = NULL;
    acc->len = 0u;
}
