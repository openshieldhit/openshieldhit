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
