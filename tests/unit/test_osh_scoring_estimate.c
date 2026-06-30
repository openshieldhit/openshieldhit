/*
 * Unit tests for osh_scoring_estimate_memory(): the allocation-free scoring
 * memory estimate that the OOM gate relies on.  Loads a committed detect.dat
 * fixture whose accumulator byte total is known exactly and checks the estimate
 * against it (including the 2x "data2" doubling for LET-average quantities).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_app_osh.h"
#include "openshieldhit/scoring.h"
#include "openshieldhit/status.h"

#ifndef OSH_TEST_FIXTURES_DIR
#define OSH_TEST_FIXTURES_DIR "."
#endif

#include "test_assert.h"

static void test_estimate_known_fixture(void) {
    char path[1024];
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_mem_estimate est;

    (void) snprintf(path, sizeof(path), "%s/scoring_estimate/detect.dat", OSH_TEST_FIXTURES_DIR);

    ASSERT_TRUE(osh_scoring_setup_from_path(path, NULL, &ws) == OSH_OK);
    ASSERT_TRUE(ws != NULL);

    ASSERT_TRUE(osh_scoring_estimate_memory(ws, &est) == OSH_OK);

    /* 100 bins: Dose = 100*8*1 = 800, DLET = 100*8*2 = 1600 -> 2400 total. */
    ASSERT_TRUE(est.npages == 2u);
    ASSERT_TRUE(est.accum_bytes == 2400u);
    ASSERT_TRUE(est.largest_page_bytes == 1600u);
    ASSERT_TRUE(strcmp(est.largest_geometry, "G") == 0);

    /* A mid-run snapshot copies only the `data` array of pages whose postprocess
     * writes data.  Dose (MeV/g) is a no-op, so it contributes nothing; DLET's
     * single data array is 100*8 = 800 (its data2 weight array is aliased, not
     * copied).  shadow_bytes is therefore strictly less than accum_bytes here. */
    ASSERT_TRUE(est.shadow_bytes == 800u);

    osh_scoring_workspace_free(ws);
}

/*
 * A degenerate geometry (an axis with 0 bins) gives every page len==0.  The
 * allocator and the shadow both round that up to one element, so the estimate
 * must too — counting a flat 0 bytes would understate what is really allocated
 * and break the "cannot drift from the real allocation" invariant the header
 * documents.  Mirrors osh_scoring_accumulator_alloc()'s `len ? len : 1`.
 */
static void test_estimate_zero_bins(void) {
    char path[1024];
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_mem_estimate est;

    (void) snprintf(path, sizeof(path), "%s/scoring_estimate/detect_zero_bins.dat", OSH_TEST_FIXTURES_DIR);

    ASSERT_TRUE(osh_scoring_setup_from_path(path, NULL, &ws) == OSH_OK);
    ASSERT_TRUE(ws != NULL);

    ASSERT_TRUE(osh_scoring_estimate_memory(ws, &est) == OSH_OK);

    /* len==0 rounds to one element: Dose = 1*8*1 = 8, DLET = 1*8*2 = 16. */
    ASSERT_TRUE(est.npages == 2u);
    ASSERT_TRUE(est.accum_bytes == 24u);
    ASSERT_TRUE(est.largest_page_bytes == 16u);
    ASSERT_TRUE(strcmp(est.largest_geometry, "G") == 0);

    /* Shadow copies DLET's single data array (rounded to one double); Dose is a
     * no-op postprocess and contributes nothing. */
    ASSERT_TRUE(est.shadow_bytes == 8u);

    osh_scoring_workspace_free(ws);
}

static void test_estimate_null_args(void) {
    struct osh_scoring_mem_estimate est;
    ASSERT_TRUE(osh_scoring_estimate_memory(NULL, &est) == OSH_EINVAL);
}

int main(void) {
    test_estimate_known_fixture();
    test_estimate_zero_bins();
    test_estimate_null_args();
    return 0;
}
