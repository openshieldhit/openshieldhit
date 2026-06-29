/*
 * Unit tests for the non-destructive snapshot shadow (issue #191).
 *
 * Coverage:
 *   test_shadow_alias_and_ownership — view aliases live; only transformed pages
 *                                     own scratch; shadow_bytes counts data-only
 *   test_snapshot_nondestructive    — snapshot leaves every live array byte-identical
 *                                     while the view holds the postprocessed values
 *   test_snapshot_reuse             — scratch is allocated once and reused across dumps
 *   test_snapshot_selector          — the output selector (G2) is passed to the sink
 *   test_snapshot_errors            — NULL sink / sink->save / shadow are rejected
 *
 * Uses a mock sink so no files are written; pure arithmetic, identical on every OS.
 */

#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "openshieldhit/const.h"
#include "openshieldhit/status.h"
#include "scoring/runtime/osh_scoring_accumulator.h"
#include "scoring/runtime/osh_scoring_runtime.h"
#include "scoring/runtime/osh_scoring_shadow.h"
#include "scoring/runtime/osh_scoring_snapshot.h"
#include "scoring/save/osh_scoring_sink.h"
#include "test_assert.h"

/* ---- A live runtime: DOSEGY (transform) + DLET (transform) + ENERGY (alias) -- */

struct live_fixture {
    struct osh_scoring_page_runtime pages[3];
    struct osh_scoring_runtime rt;
};

static void
make_page(struct osh_scoring_page_runtime *p, enum osh_scoring_score_kind kind, size_t len, int want_data2) {
    memset(p, 0, sizeof(*p));
    ASSERT_TRUE(osh_scoring_accumulator_alloc(&p->acc, len, want_data2) == OSH_OK);
    p->len = len;
    p->score_kind = kind;
    p->has_data2 = (char) (want_data2 ? 1 : 0);
    p->postproc = (kind == OSH_SCORING_SCORE_DLET) ? OSH_SCORING_POSTPROC_AVER : OSH_SCORING_POSTPROC_NONE;
}

static void live_fixture_init(struct live_fixture *f) {
    make_page(&f->pages[0], OSH_SCORING_SCORE_DOSEGY, 3u, 0);
    f->pages[0].acc.data[0] = 1.0;
    f->pages[0].acc.data[1] = 2.0;
    f->pages[0].acc.data[2] = 3.0;

    make_page(&f->pages[1], OSH_SCORING_SCORE_DLET, 3u, 1);
    f->pages[1].acc.data[0] = 10.0;
    f->pages[1].acc.data2[0] = 2.0; /* -> 5 */
    f->pages[1].acc.data[1] = 9.0;
    f->pages[1].acc.data2[1] = 3.0; /* -> 3 */
    f->pages[1].acc.data[2] = 7.0;
    f->pages[1].acc.data2[2] = 0.0; /* weight ~0 -> 0 */

    make_page(&f->pages[2], OSH_SCORING_SCORE_ENERGY, 3u, 0);
    f->pages[2].acc.data[0] = 100.0;
    f->pages[2].acc.data[1] = 200.0;
    f->pages[2].acc.data[2] = 300.0;

    memset(&f->rt, 0, sizeof(f->rt));
    f->rt.pages = f->pages;
    f->rt.npages = 3u;
}

static void live_fixture_free(struct live_fixture *f) {
    size_t i;
    for (i = 0; i < 3u; ++i) {
        osh_scoring_accumulator_free(&f->pages[i].acc);
    }
}

/* ---- Mock sink ------------------------------------------------------------ */

struct mock_sink_ctx {
    int calls;
    struct osh_scoring_runtime const *last_rt;
    unsigned long long last_nstat;
    size_t const *last_want;
    size_t last_n_want;
};

static enum osh_status mock_save(
    void *ctx, struct osh_scoring_runtime const *rt, unsigned long long nstat, size_t const *want, size_t n_want) {
    struct mock_sink_ctx *m = (struct mock_sink_ctx *) ctx;
    m->calls += 1;
    m->last_rt = rt;
    m->last_nstat = nstat;
    m->last_want = want;
    m->last_n_want = n_want;
    return OSH_OK;
}

/* ---- view aliases live; only transformed pages own scratch ---------------- */

static void test_shadow_alias_and_ownership(void) {
    struct live_fixture f;
    struct osh_scoring_shadow shadow;

    live_fixture_init(&f);
    ASSERT_TRUE(osh_scoring_shadow_init(&shadow, &f.rt) == OSH_OK);

    /* Distinct page array, but every page's data still aliases live before refresh. */
    ASSERT_TRUE(shadow.view.pages != f.rt.pages);
    ASSERT_TRUE(shadow.view.npages == 3u);
    ASSERT_TRUE(shadow.pages[0].acc.data == f.pages[0].acc.data);
    ASSERT_TRUE(shadow.pages[1].acc.data == f.pages[1].acc.data);
    ASSERT_TRUE(shadow.pages[2].acc.data == f.pages[2].acc.data);

    /* data-only over transformed pages: DOSEGY(3) + DLET(3) = 6 doubles = 48 B. */
    ASSERT_TRUE(osh_scoring_shadow_bytes(&shadow) == 48u);

    ASSERT_TRUE(osh_scoring_shadow_refresh(&shadow) == OSH_OK);

    /* After refresh: transformed pages own a private buffer; the simple page still
     * aliases live, and data2 is aliased on every page (never copied). */
    ASSERT_TRUE(shadow.pages[0].acc.data != f.pages[0].acc.data);
    ASSERT_TRUE(shadow.pages[1].acc.data != f.pages[1].acc.data);
    ASSERT_TRUE(shadow.pages[2].acc.data == f.pages[2].acc.data);
    ASSERT_TRUE(shadow.pages[1].acc.data2 == f.pages[1].acc.data2);
    ASSERT_TRUE(shadow.scratch[0] != NULL);
    ASSERT_TRUE(shadow.scratch[1] != NULL);
    ASSERT_TRUE(shadow.scratch[2] == NULL);

    osh_scoring_shadow_free(&shadow);
    live_fixture_free(&f);
}

/* ---- Headline: snapshot is non-destructive -------------------------------- */

static void test_snapshot_nondestructive(void) {
    struct live_fixture f;
    struct osh_scoring_shadow shadow;
    struct mock_sink_ctx ctx = {0};
    struct osh_scoring_sink sink;
    struct osh_scoring_runtime const *view;

    /* Byte images of live state before the snapshot. */
    double d0[3], d1[3], w1[3], d2[3];

    live_fixture_init(&f);
    memcpy(d0, f.pages[0].acc.data, sizeof(d0));
    memcpy(d1, f.pages[1].acc.data, sizeof(d1));
    memcpy(w1, f.pages[1].acc.data2, sizeof(w1));
    memcpy(d2, f.pages[2].acc.data, sizeof(d2));

    sink.save = mock_save;
    sink.ctx = &ctx;

    ASSERT_TRUE(osh_scoring_shadow_init(&shadow, &f.rt) == OSH_OK);
    ASSERT_TRUE(osh_scoring_snapshot_save(&sink, &shadow, 12345ull, NULL, 0u) == OSH_OK);

    /* The sink saw the shadow view and the completed count. */
    ASSERT_TRUE(ctx.calls == 1);
    ASSERT_TRUE(ctx.last_rt == osh_scoring_shadow_view(&shadow));
    ASSERT_TRUE(ctx.last_nstat == 12345ull);

    /* View holds the postprocessed values. */
    view = ctx.last_rt;
    ASSERT_TRUE(view->pages[0].acc.data[0] == 1.0 * OSH_MEVG2GY);
    ASSERT_TRUE(view->pages[0].acc.data[2] == 3.0 * OSH_MEVG2GY);
    ASSERT_TRUE(view->pages[1].acc.data[0] == 5.0);
    ASSERT_TRUE(view->pages[1].acc.data[1] == 3.0);
    ASSERT_TRUE(view->pages[1].acc.data[2] == 0.0);
    ASSERT_TRUE(view->pages[1].has_data2 == 0);
    ASSERT_TRUE(view->pages[2].acc.data[1] == 200.0); /* simple page passes raw values through */

    /* Live state is byte-identical — the whole point. */
    ASSERT_TRUE(memcmp(f.pages[0].acc.data, d0, sizeof(d0)) == 0);
    ASSERT_TRUE(memcmp(f.pages[1].acc.data, d1, sizeof(d1)) == 0);
    ASSERT_TRUE(memcmp(f.pages[1].acc.data2, w1, sizeof(w1)) == 0);
    ASSERT_TRUE(memcmp(f.pages[2].acc.data, d2, sizeof(d2)) == 0);
    ASSERT_TRUE(f.pages[1].has_data2 == 1); /* live flag untouched */

    osh_scoring_shadow_free(&shadow);
    live_fixture_free(&f);
}

/* ---- Scratch allocated once, reused across dumps -------------------------- */

static void test_snapshot_reuse(void) {
    struct live_fixture f;
    struct osh_scoring_shadow shadow;
    struct mock_sink_ctx ctx = {0};
    struct osh_scoring_sink sink;
    double *s0, *s1;
    double d1[3], w1[3];

    live_fixture_init(&f);
    memcpy(d1, f.pages[1].acc.data, sizeof(d1));
    memcpy(w1, f.pages[1].acc.data2, sizeof(w1));

    sink.save = mock_save;
    sink.ctx = &ctx;

    ASSERT_TRUE(osh_scoring_shadow_init(&shadow, &f.rt) == OSH_OK);
    ASSERT_TRUE(osh_scoring_snapshot_save(&sink, &shadow, 10ull, NULL, 0u) == OSH_OK);
    s0 = shadow.scratch[0];
    s1 = shadow.scratch[1];

    /* Second dump must not reallocate scratch. */
    ASSERT_TRUE(osh_scoring_snapshot_save(&sink, &shadow, 20ull, NULL, 0u) == OSH_OK);
    ASSERT_TRUE(shadow.scratch[0] == s0);
    ASSERT_TRUE(shadow.scratch[1] == s1);
    ASSERT_TRUE(ctx.calls == 2);
    ASSERT_TRUE(ctx.last_nstat == 20ull);

    /* Live still intact after two dumps. */
    ASSERT_TRUE(memcmp(f.pages[1].acc.data, d1, sizeof(d1)) == 0);
    ASSERT_TRUE(memcmp(f.pages[1].acc.data2, w1, sizeof(w1)) == 0);

    osh_scoring_shadow_free(&shadow);
    live_fixture_free(&f);
}

/* ---- Output selector (G2) is passed through to the sink ------------------- */

static void test_snapshot_selector(void) {
    struct live_fixture f;
    struct osh_scoring_shadow shadow;
    struct mock_sink_ctx ctx = {0};
    struct osh_scoring_sink sink;
    size_t const want[1] = {2u};

    live_fixture_init(&f);
    sink.save = mock_save;
    sink.ctx = &ctx;

    ASSERT_TRUE(osh_scoring_shadow_init(&shadow, &f.rt) == OSH_OK);
    ASSERT_TRUE(osh_scoring_snapshot_save(&sink, &shadow, 7ull, want, 1u) == OSH_OK);
    ASSERT_TRUE(ctx.last_want == want);
    ASSERT_TRUE(ctx.last_n_want == 1u);

    osh_scoring_shadow_free(&shadow);
    live_fixture_free(&f);
}

/* ---- Error paths ---------------------------------------------------------- */

static void test_snapshot_errors(void) {
    struct live_fixture f;
    struct osh_scoring_shadow shadow;
    struct mock_sink_ctx ctx = {0};
    struct osh_scoring_sink good;
    struct osh_scoring_sink no_save;

    live_fixture_init(&f);
    good.save = mock_save;
    good.ctx = &ctx;
    no_save.save = NULL;
    no_save.ctx = &ctx;

    ASSERT_TRUE(osh_scoring_shadow_init(&shadow, &f.rt) == OSH_OK);
    ASSERT_TRUE(osh_scoring_snapshot_save(NULL, &shadow, 1ull, NULL, 0u) == OSH_EINVAL);
    ASSERT_TRUE(osh_scoring_snapshot_save(&no_save, &shadow, 1ull, NULL, 0u) == OSH_EINVAL);
    ASSERT_TRUE(osh_scoring_snapshot_save(&good, NULL, 1ull, NULL, 0u) == OSH_EINVAL);
    ASSERT_TRUE(ctx.calls == 0); /* nothing was saved */

    osh_scoring_shadow_free(&shadow);
    live_fixture_free(&f);
}

int main(void) {
    test_shadow_alias_and_ownership();
    test_snapshot_nondestructive();
    test_snapshot_reuse();
    test_snapshot_selector();
    test_snapshot_errors();
    printf("All osh_scoring_shadow tests passed.\n");
    return 0;
}
