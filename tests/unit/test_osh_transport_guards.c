#include <string.h>

#include "test_assert.h"
#include "transport/osh_transport.h"
#include "transport/osh_transport_ion.h"

/*
 * Argument-validation guards for the transport entry points.  These exercise the
 * early OSH_EINVAL returns without building a full runtime: the borrowed runtime
 * pointers are non-NULL sentinels that the guards reject on before ever
 * dereferencing them, so a valid, never-read address is all that is required.
 */

static void test_run_minimal_rejects_bad_args(void) {
    struct osh_transport_context ctx;

    /* NULL context. */
    ASSERT_TRUE(osh_transport_run_minimal(NULL, NULL, NULL, NULL, NULL) == OSH_EINVAL);

    /* nstat == 0 is rejected before any pool is touched, so a zeroed context
     * (NULL pools) is enough to reach and return from that guard. */
    memset(&ctx, 0, sizeof(ctx));
    ctx.params.nstat = 0u;
    ASSERT_TRUE(osh_transport_run_minimal(&ctx, NULL, NULL, NULL, NULL) == OSH_EINVAL);
}

static void test_ion_run_range_validates_arguments(void) {
    struct osh_transport_context ctx;
    double storage = 0.0;
    void *nn = &storage; /* non-NULL sentinel; the guards never dereference it */
    size_t completed = 12345u;

    /* NULL context: rejected, and completed_out is zeroed up front. */
    ASSERT_TRUE(osh_transport_ion_run_range(NULL, nn, nn, nn, nn, 0u, 10u, &completed) == OSH_EINVAL);
    ASSERT_TRUE(completed == 0u);

    /* Missing pool / scratch on the context. */
    memset(&ctx, 0, sizeof(ctx));
    ASSERT_TRUE(osh_transport_ion_run_range(&ctx, nn, nn, nn, nn, 0u, 10u, &completed) == OSH_EINVAL);

    /* Wire non-NULL sentinel pool + scratch so the later guards are reachable. */
    ctx.ion_pool = nn;
    ctx.zone_refs = nn;
    ctx.dist_batch = nn;

    /* Empty / inverted range (hist_hi <= hist_lo). */
    ASSERT_TRUE(osh_transport_ion_run_range(&ctx, nn, nn, nn, nn, 5u, 5u, &completed) == OSH_EINVAL);

    /* nstat == 0. */
    ctx.params.nstat = 0u;
    ASSERT_TRUE(osh_transport_ion_run_range(&ctx, nn, nn, nn, nn, 0u, 10u, &completed) == OSH_EINVAL);

    /* Out-of-range DELTAE. */
    ctx.params.nstat = 10u;
    ctx.params.deltae = 0.0f;
    ASSERT_TRUE(osh_transport_ion_run_range(&ctx, nn, nn, nn, nn, 0u, 10u, &completed) == OSH_EINVAL);
}

int main(void) {
    test_run_minimal_rejects_bad_args();
    test_ion_run_range_validates_arguments();
    return 0;
}
