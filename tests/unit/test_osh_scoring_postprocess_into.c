/*
 * Unit tests for the out-of-place postprocess (issue #191).
 *
 * Coverage:
 *   test_writes_data_predicate    — only DOSEGY + LET/Qeff write data
 *   test_inplace_dosegy           — in-place wrapper rescales MeV/g -> Gy
 *   test_inplace_let              — in-place wrapper finalises data/data2, clears flags
 *   test_inplace_simple_noop      — simple scorers are untouched
 *   test_into_out_of_place        — dst gets the transform, src stays unchanged
 *   test_into_let_nondestructive  — out-of-place LET leaves src data + data2 intact
 *   test_into_errors              — NULL args and npages/len mismatch are rejected
 *
 * Pure arithmetic over hand-built runtimes, so identical on every OS.
 */

#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "openshieldhit/const.h"
#include "openshieldhit/status.h"
#include "scoring/runtime/osh_scoring_accumulator.h"
#include "scoring/runtime/osh_scoring_postprocess.h"
#include "scoring/runtime/osh_scoring_runtime.h"
#include "test_assert.h"

/* Build a single page descriptor with an allocated accumulator. */
static void
make_page(struct osh_scoring_page_runtime *p, enum osh_scoring_score_kind kind, size_t len, int want_data2) {
    memset(p, 0, sizeof(*p));
    ASSERT_TRUE(osh_scoring_accumulator_alloc(&p->acc, len, want_data2) == OSH_OK);
    p->len = len;
    p->score_kind = kind;
    p->has_data2 = (char) (want_data2 ? 1 : 0);
    p->divide = 0;
    p->postproc = (kind == OSH_SCORING_SCORE_DLET || kind == OSH_SCORING_SCORE_TLET || kind == OSH_SCORING_SCORE_DQEFF
                   || kind == OSH_SCORING_SCORE_TQEFF)
                      ? OSH_SCORING_POSTPROC_AVER
                      : OSH_SCORING_POSTPROC_NONE;
}

static void make_runtime(struct osh_scoring_runtime *rt, struct osh_scoring_page_runtime *pages, size_t npages) {
    memset(rt, 0, sizeof(*rt));
    rt->pages = pages;
    rt->npages = npages;
}

/* ---- The predicate that decides which pages need a private buffer --------- */

static void test_writes_data_predicate(void) {
    ASSERT_TRUE(osh_scoring_postprocess_writes_data(OSH_SCORING_SCORE_DOSEGY) == 1);
    ASSERT_TRUE(osh_scoring_postprocess_writes_data(OSH_SCORING_SCORE_DLET) == 1);
    ASSERT_TRUE(osh_scoring_postprocess_writes_data(OSH_SCORING_SCORE_TLET) == 1);
    ASSERT_TRUE(osh_scoring_postprocess_writes_data(OSH_SCORING_SCORE_DQEFF) == 1);
    ASSERT_TRUE(osh_scoring_postprocess_writes_data(OSH_SCORING_SCORE_TQEFF) == 1);

    /* DOSE/FLUENCE (and NKERMA) now transform in postprocess too — they divide by
     * the per-bin volume, so they write data out-of-place like DOSEGY/LET. */
    ASSERT_TRUE(osh_scoring_postprocess_writes_data(OSH_SCORING_SCORE_DOSE) == 1);
    ASSERT_TRUE(osh_scoring_postprocess_writes_data(OSH_SCORING_SCORE_FLUENCE) == 1);
    ASSERT_TRUE(osh_scoring_postprocess_writes_data(OSH_SCORING_SCORE_NKERMA) == 1);

    ASSERT_TRUE(osh_scoring_postprocess_writes_data(OSH_SCORING_SCORE_ENERGY) == 0);
    ASSERT_TRUE(osh_scoring_postprocess_writes_data(OSH_SCORING_SCORE_COUNT) == 0);
    ASSERT_TRUE(osh_scoring_postprocess_writes_data(OSH_SCORING_SCORE_UNKNOWN) == 0);
}

/* ---- In-place wrapper: DOSEGY rescale ------------------------------------- */

static void test_inplace_dosegy(void) {
    struct osh_scoring_page_runtime page;
    struct osh_scoring_runtime rt;
    size_t i;

    make_page(&page, OSH_SCORING_SCORE_DOSEGY, 3u, 0);
    page.acc.data[0] = 1.0;
    page.acc.data[1] = 2.0;
    page.acc.data[2] = 3.0;
    make_runtime(&rt, &page, 1u);

    ASSERT_TRUE(osh_scoring_postprocess(&rt) == OSH_OK);
    for (i = 0; i < 3u; ++i) {
        ASSERT_TRUE(page.acc.data[i] == (double) (i + 1) * OSH_MEVG2GY);
    }

    osh_scoring_accumulator_free(&page.acc);
}

/* ---- In-place wrapper: LET finalise + flag clearing ----------------------- */

static void test_inplace_let(void) {
    struct osh_scoring_page_runtime page;
    struct osh_scoring_runtime rt;

    make_page(&page, OSH_SCORING_SCORE_DLET, 3u, 1);
    page.acc.data[0] = 10.0;
    page.acc.data2[0] = 2.0; /* -> 5 */
    page.acc.data[1] = 9.0;
    page.acc.data2[1] = 3.0; /* -> 3 */
    page.acc.data[2] = 7.0;
    page.acc.data2[2] = 0.0; /* weight ~0 -> 0 */
    make_runtime(&rt, &page, 1u);

    ASSERT_TRUE(osh_scoring_postprocess(&rt) == OSH_OK);
    ASSERT_TRUE(page.acc.data[0] == 5.0);
    ASSERT_TRUE(page.acc.data[1] == 3.0);
    ASSERT_TRUE(page.acc.data[2] == 0.0);
    ASSERT_TRUE(page.has_data2 == 0);
    ASSERT_TRUE(page.divide == 0);

    osh_scoring_accumulator_free(&page.acc);
}

/* ---- In-place wrapper: simple scorer is a no-op --------------------------- */

static void test_inplace_simple_noop(void) {
    struct osh_scoring_page_runtime page;
    struct osh_scoring_runtime rt;

    make_page(&page, OSH_SCORING_SCORE_ENERGY, 3u, 0);
    page.acc.data[0] = 11.0;
    page.acc.data[1] = 22.0;
    page.acc.data[2] = 33.0;
    make_runtime(&rt, &page, 1u);

    ASSERT_TRUE(osh_scoring_postprocess(&rt) == OSH_OK);
    ASSERT_TRUE(page.acc.data[0] == 11.0);
    ASSERT_TRUE(page.acc.data[1] == 22.0);
    ASSERT_TRUE(page.acc.data[2] == 33.0);

    osh_scoring_accumulator_free(&page.acc);
}

/* ---- In-place wrapper: single-shot guard --------------------------------- */

static void test_inplace_single_shot_guard(void) {
    struct osh_scoring_page_runtime page;
    struct osh_scoring_runtime rt;
    size_t i;

    make_page(&page, OSH_SCORING_SCORE_DOSEGY, 3u, 0);
    page.acc.data[0] = 1.0;
    page.acc.data[1] = 2.0;
    page.acc.data[2] = 3.0;
    make_runtime(&rt, &page, 1u);

    /* First call finalises in place. */
    ASSERT_TRUE(osh_scoring_postprocess(&rt) == OSH_OK);
    for (i = 0; i < 3u; ++i) {
        ASSERT_TRUE(page.acc.data[i] == (double) (i + 1) * OSH_MEVG2GY);
    }

    /* A second in-place call is refused (would double-apply the MeV/g->Gy
     * rescale); the buffers stay at their once-finalised values. */
    ASSERT_TRUE(osh_scoring_postprocess(&rt) == OSH_ESTATE);
    for (i = 0; i < 3u; ++i) {
        ASSERT_TRUE(page.acc.data[i] == (double) (i + 1) * OSH_MEVG2GY);
    }

    osh_scoring_accumulator_free(&page.acc);
}

/* ---- Direct in-place _into(rt, rt): same single-shot guard as the wrapper -- */

static void test_into_inplace_single_shot_guard(void) {
    struct osh_scoring_page_runtime page;
    struct osh_scoring_runtime rt;
    size_t i;

    make_page(&page, OSH_SCORING_SCORE_DOSEGY, 3u, 0);
    page.acc.data[0] = 1.0;
    page.acc.data[1] = 2.0;
    page.acc.data[2] = 3.0;
    make_runtime(&rt, &page, 1u);

    /* Calling the primitive in place (dst == src) must be guarded exactly like the
     * osh_scoring_postprocess() wrapper: the first call finalises, the second is
     * refused with OSH_ESTATE rather than double-applying the MeV/g->Gy rescale. */
    ASSERT_TRUE(osh_scoring_postprocess_into(&rt, &rt) == OSH_OK);
    for (i = 0; i < 3u; ++i) {
        ASSERT_TRUE(page.acc.data[i] == (double) (i + 1) * OSH_MEVG2GY);
    }
    ASSERT_TRUE(rt.postprocessed == 1);

    ASSERT_TRUE(osh_scoring_postprocess_into(&rt, &rt) == OSH_ESTATE);
    for (i = 0; i < 3u; ++i) {
        ASSERT_TRUE(page.acc.data[i] == (double) (i + 1) * OSH_MEVG2GY);
    }

    /* An out-of-place call after an in-place finalisation is still allowed — the
     * guard is in-place-only, so the shadow/dump path stays repeatable. */
    {
        struct osh_scoring_page_runtime dpage;
        struct osh_scoring_runtime dst;

        make_page(&dpage, OSH_SCORING_SCORE_DOSEGY, 3u, 0);
        make_runtime(&dst, &dpage, 1u);
        ASSERT_TRUE(osh_scoring_postprocess_into(&dst, &rt) == OSH_OK);
        ASSERT_TRUE(dst.postprocessed == 0); /* out-of-place never sets the flag */
        osh_scoring_accumulator_free(&dpage.acc);
    }

    osh_scoring_accumulator_free(&page.acc);
}

/* ---- Out-of-place: dst transformed, src untouched ------------------------- */

static void test_into_out_of_place(void) {
    struct osh_scoring_page_runtime spage;
    struct osh_scoring_page_runtime dpage;
    struct osh_scoring_runtime src;
    struct osh_scoring_runtime dst;
    size_t i;

    make_page(&spage, OSH_SCORING_SCORE_DOSEGY, 3u, 0);
    make_page(&dpage, OSH_SCORING_SCORE_DOSEGY, 3u, 0);
    for (i = 0; i < 3u; ++i) {
        spage.acc.data[i] = (double) (i + 1);
        dpage.acc.data[i] = -1.0; /* sentinel: must be overwritten */
    }
    make_runtime(&src, &spage, 1u);
    make_runtime(&dst, &dpage, 1u);

    ASSERT_TRUE(osh_scoring_postprocess_into(&dst, &src) == OSH_OK);
    for (i = 0; i < 3u; ++i) {
        ASSERT_TRUE(dpage.acc.data[i] == (double) (i + 1) * OSH_MEVG2GY); /* dst transformed */
        ASSERT_TRUE(spage.acc.data[i] == (double) (i + 1));               /* src untouched */
    }

    osh_scoring_accumulator_free(&spage.acc);
    osh_scoring_accumulator_free(&dpage.acc);
}

/* ---- Out-of-place LET: src data AND data2 left intact --------------------- */

static void test_into_let_nondestructive(void) {
    struct osh_scoring_page_runtime spage;
    struct osh_scoring_page_runtime dpage;
    struct osh_scoring_runtime src;
    struct osh_scoring_runtime dst;

    make_page(&spage, OSH_SCORING_SCORE_DLET, 2u, 1);
    make_page(&dpage, OSH_SCORING_SCORE_DLET, 2u, 1);
    spage.acc.data[0] = 12.0;
    spage.acc.data2[0] = 3.0; /* -> 4 */
    spage.acc.data[1] = 10.0;
    spage.acc.data2[1] = 5.0; /* -> 2 */
    make_runtime(&src, &spage, 1u);
    make_runtime(&dst, &dpage, 1u);

    ASSERT_TRUE(osh_scoring_postprocess_into(&dst, &src) == OSH_OK);

    /* dst holds the finalised average with flags cleared. */
    ASSERT_TRUE(dpage.acc.data[0] == 4.0);
    ASSERT_TRUE(dpage.acc.data[1] == 2.0);
    ASSERT_TRUE(dpage.has_data2 == 0);
    ASSERT_TRUE(dpage.divide == 0);

    /* src is byte-identical: data is the raw weighted sum, data2 the weight, and
     * its flags are unchanged so live accumulation could continue. */
    ASSERT_TRUE(spage.acc.data[0] == 12.0);
    ASSERT_TRUE(spage.acc.data[1] == 10.0);
    ASSERT_TRUE(spage.acc.data2[0] == 3.0);
    ASSERT_TRUE(spage.acc.data2[1] == 5.0);
    ASSERT_TRUE(spage.has_data2 == 1);

    osh_scoring_accumulator_free(&spage.acc);
    osh_scoring_accumulator_free(&dpage.acc);
}

/* ---- Error paths ---------------------------------------------------------- */

static void test_into_errors(void) {
    struct osh_scoring_page_runtime pa;
    struct osh_scoring_page_runtime pb;
    struct osh_scoring_runtime a;
    struct osh_scoring_runtime b;
    struct osh_scoring_runtime empty;

    make_page(&pa, OSH_SCORING_SCORE_ENERGY, 3u, 0);
    make_page(&pb, OSH_SCORING_SCORE_ENERGY, 4u, 0); /* different len */
    make_runtime(&a, &pa, 1u);
    make_runtime(&b, &pb, 1u);
    make_runtime(&empty, NULL, 0u);

    ASSERT_TRUE(osh_scoring_postprocess_into(NULL, &a) == OSH_EINVAL);
    ASSERT_TRUE(osh_scoring_postprocess_into(&a, NULL) == OSH_EINVAL);
    ASSERT_TRUE(osh_scoring_postprocess_into(&a, &empty) == OSH_EINVAL); /* npages mismatch */
    ASSERT_TRUE(osh_scoring_postprocess_into(&a, &b) == OSH_EINVAL);     /* per-page len mismatch */

    osh_scoring_accumulator_free(&pa.acc);
    osh_scoring_accumulator_free(&pb.acc);

    /* npages > 0 with a NULL page array is rejected up front. */
    {
        struct osh_scoring_runtime bad_pages;
        struct osh_scoring_runtime ok;
        struct osh_scoring_page_runtime pc;

        make_page(&pc, OSH_SCORING_SCORE_ENERGY, 2u, 0);
        make_runtime(&ok, &pc, 1u);
        memset(&bad_pages, 0, sizeof(bad_pages));
        bad_pages.npages = 1u; /* pages stays NULL */

        ASSERT_TRUE(osh_scoring_postprocess_into(&bad_pages, &ok) == OSH_EINVAL);
        osh_scoring_accumulator_free(&pc.acc);
    }
}

/* dst/src describing different kinds is rejected (mispaired pages). */
static void test_into_kind_mismatch(void) {
    struct osh_scoring_page_runtime spage;
    struct osh_scoring_page_runtime dpage;
    struct osh_scoring_runtime src;
    struct osh_scoring_runtime dst;

    make_page(&spage, OSH_SCORING_SCORE_DOSEGY, 3u, 0);
    make_page(&dpage, OSH_SCORING_SCORE_DLET, 3u, 1); /* same len, different kind */
    make_runtime(&src, &spage, 1u);
    make_runtime(&dst, &dpage, 1u);

    ASSERT_TRUE(osh_scoring_postprocess_into(&dst, &src) == OSH_EINVAL);

    osh_scoring_accumulator_free(&spage.acc);
    osh_scoring_accumulator_free(&dpage.acc);
}

/* Missing required arrays are rejected, not dereferenced. */
static void test_into_missing_arrays(void) {
    /* LET with no data2 (the divisor) -> EINVAL. */
    {
        struct osh_scoring_page_runtime page;
        struct osh_scoring_runtime rt;
        make_page(&page, OSH_SCORING_SCORE_DLET, 3u, 0); /* want_data2 = 0 -> data2 NULL */
        make_runtime(&rt, &page, 1u);
        ASSERT_TRUE(osh_scoring_postprocess(&rt) == OSH_EINVAL);
        osh_scoring_accumulator_free(&page.acc);
    }
    /* DOSEGY with a NULL dst data buffer -> EINVAL. */
    {
        struct osh_scoring_page_runtime spage;
        struct osh_scoring_page_runtime dpage;
        struct osh_scoring_runtime src;
        struct osh_scoring_runtime dst;

        make_page(&spage, OSH_SCORING_SCORE_DOSEGY, 3u, 0);
        memset(&dpage, 0, sizeof(dpage)); /* acc.data NULL */
        dpage.acc.len = 3u;
        dpage.len = 3u;
        dpage.score_kind = OSH_SCORING_SCORE_DOSEGY;
        make_runtime(&src, &spage, 1u);
        make_runtime(&dst, &dpage, 1u);

        ASSERT_TRUE(osh_scoring_postprocess_into(&dst, &src) == OSH_EINVAL);
        osh_scoring_accumulator_free(&spage.acc);
    }
}

/* Unhandled postproc shapes return OSH_ENOTSUP. */
static void test_into_unsupported(void) {
    /* A simple scorer carrying has_data2 is an unhandled two-pass page. */
    {
        struct osh_scoring_page_runtime page;
        struct osh_scoring_runtime rt;
        make_page(&page, OSH_SCORING_SCORE_ENERGY, 2u, 0);
        page.has_data2 = 1;
        make_runtime(&rt, &page, 1u);
        ASSERT_TRUE(osh_scoring_postprocess(&rt) == OSH_ENOTSUP);
        osh_scoring_accumulator_free(&page.acc);
    }
    /* AVER postproc on a non-LET kind has no finaliser. */
    {
        struct osh_scoring_page_runtime page;
        struct osh_scoring_runtime rt;
        make_page(&page, OSH_SCORING_SCORE_ENERGY, 2u, 0);
        page.postproc = OSH_SCORING_POSTPROC_AVER;
        make_runtime(&rt, &page, 1u);
        ASSERT_TRUE(osh_scoring_postprocess(&rt) == OSH_ENOTSUP);
        osh_scoring_accumulator_free(&page.acc);
    }
}

int main(void) {
    test_writes_data_predicate();
    test_inplace_dosegy();
    test_inplace_let();
    test_inplace_simple_noop();
    test_inplace_single_shot_guard();
    test_into_inplace_single_shot_guard();
    test_into_out_of_place();
    test_into_let_nondestructive();
    test_into_errors();
    test_into_kind_mismatch();
    test_into_missing_arrays();
    test_into_unsupported();
    printf("All osh_scoring_postprocess_into tests passed.\n");
    return 0;
}
