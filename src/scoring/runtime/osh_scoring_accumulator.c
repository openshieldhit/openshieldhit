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
    acc->weight = 0.0;
    acc->nbatch = 0u;

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
    if (!acc) {
        return;
    }
    /* Reset the batch bookkeeping unconditionally: alloc() supports len==0 (it
     * still allocates one element), so a zeroed accumulator must read as "no
     * batches yet" even then. */
    acc->weight = 0.0;
    acc->nbatch = 0u;
    if (acc->len == 0u) {
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

/*
 * Numerically-stable parallel merge of two Welford M2 arrays, in place: dst_m2
 * becomes the combined M2 of batch A (dst_*, holding running sum dst_sum and
 * weight dst_w) and batch B (src_*).  Per-bin means are derived as sum/weight;
 * the weights are the batches' history counts.
 *
 *   M2 = M2_A + M2_B + (mean_B − mean_A)² · w_A·w_B / (w_A + w_B)
 *
 * The cross-term is what a plain element-wise += cannot express; dropping it
 * silently corrupts the variance.  This must run before the additive fields
 * (dst_sum / dst_w) are folded, since it reads dst's pre-merge values.  Empty
 * batches act as the identity.
 *
 * M2 is the running sum-of-squared-deviations of Welford's online variance
 * (Welford 1962); the weighted form is West (1979); the pairwise combine below is
 * Schubert & Gertz (2018).  Full citations in osh_scoring_accumulator.h.
 */
static void merge_m2(double *dst_m2,
                     double const *dst_sum,
                     double dst_w,
                     double const *src_m2,
                     double const *src_sum,
                     double src_w,
                     size_t len) {
    size_t i;
    double inv_dst_w;
    double inv_src_w;
    double cross;
    if (!dst_m2 || !src_m2 || !dst_sum || !src_sum) {
        return;
    }
    if (src_w == 0.0) {
        return; /* B empty: A unchanged */
    }
    if (dst_w == 0.0) {
        /* A empty: the combined batch is B. */
        for (i = 0; i < len; ++i) {
            dst_m2[i] = src_m2[i];
        }
        return;
    }
    /* All three factors are loop-invariant; hoist them so the per-bin work is
     * multiplies only (no divisions), which matters for large meshes. */
    inv_dst_w = 1.0 / dst_w;
    inv_src_w = 1.0 / src_w;
    cross = dst_w * src_w / (dst_w + src_w);
    for (i = 0; i < len; ++i) {
        double const mean_a = dst_sum[i] * inv_dst_w;
        double const mean_b = src_sum[i] * inv_src_w;
        double const delta = mean_b - mean_a;
        dst_m2[i] += src_m2[i] + delta * delta * cross;
    }
}

/* A variance-tracking accumulator must be self-consistent, or merge_m2() would
 * silently short-circuit and produce wrong M2 state: an M2 array needs its
 * companion sum array (it derives the mean from it), and a batch's weight and
 * count must agree on emptiness (W == 0 iff B == 0). */
static int variance_inconsistent(struct osh_scoring_accumulator const *a) {
    if (a->data_var && !a->data) {
        return 1;
    }
    if (a->data2_var && !a->data2) {
        return 1;
    }
    if ((a->data_var || a->data2_var) && ((a->weight == 0.0) != (a->nbatch == 0u))) {
        return 1;
    }
    return 0;
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
    /* Reject inconsistent variance bookkeeping rather than silently mis-merging. */
    if (variance_inconsistent(dst) || variance_inconsistent(src)) {
        return OSH_EINVAL;
    }
    /* Variance (Welford M2) arrays first — they need dst's pre-merge sums and
     * weight to form the Schubert-Gertz cross-term.  data_var pairs with data,
     * data2_var with data2; both use the per-accumulator history weight.  NULL
     * variance arrays (the current default) make these no-ops, leaving the plain
     * additive reduction over data / data2. */
    merge_m2(dst->data_var, dst->data, dst->weight, src->data_var, src->data, src->weight, dst->len);
    merge_m2(dst->data2_var, dst->data2, dst->weight, src->data2_var, src->data2, src->weight, dst->len);

    /* Additive fields: raw sums and the batch bookkeeping. */
    merge_array(dst->data, src->data, dst->len);
    merge_array(dst->data2, src->data2, dst->len);
    dst->weight += src->weight;
    dst->nbatch += src->nbatch;
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
    acc->weight = 0.0;
    acc->nbatch = 0u;
}
