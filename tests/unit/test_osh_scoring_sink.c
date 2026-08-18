/*
 * Unit tests for the native file sink (G1) and the output selector (G2),
 * exercised end-to-end on a real compiled workspace/runtime (issue #191).
 *
 * Coverage:
 *   test_save_outputs_errors  — NULL / nstat==0 / count-mismatch / bad-index rejects
 *   test_file_sink_selector   — file sink writes only the selected output, then all
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_app_osh.h"
#include "openshieldhit/scoring.h"
#include "openshieldhit/status.h"
#include "scoring/runtime/osh_scoring_compile.h"
#include "scoring/runtime/osh_scoring_shadow.h"
#include "scoring/runtime/osh_scoring_snapshot.h"
#include "scoring/save/osh_scoring_save.h"
#include "scoring/save/osh_scoring_sink.h"
#include "test_assert.h"

#define DETECT_PATH "osh_scoring_sink_detect.tmp"
#define OUT0 "sink_out0.txt"
#define OUT1 "sink_out1.txt"

static char const *const DETECT_TEXT = "Geometry Mesh\n"
                                       "    Name G\n"
                                       "    X 0 2 2\n"
                                       "    Y 0 1 1\n"
                                       "    Z 0 1 1\n"
                                       "\n"
                                       "Output\n"
                                       "    Filename " OUT0 "\n"
                                       "    FileFormat ASCII\n"
                                       "    Geo G\n"
                                       "    Quantity ENERGY\n"
                                       "\n"
                                       "Output\n"
                                       "    Filename " OUT1 "\n"
                                       "    FileFormat ASCII\n"
                                       "    Geo G\n"
                                       "    Quantity ENERGY\n";

static void write_detect_file(void) {
    FILE *fp = fopen(DETECT_PATH, "w");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fputs(DETECT_TEXT, fp) >= 0);
    ASSERT_TRUE(fclose(fp) == 0);
}

static int file_exists(char const *path) {
    FILE *fp = fopen(path, "r");
    if (fp) {
        (void) fclose(fp);
        return 1;
    }
    return 0;
}

/* Build a real two-output runtime with seeded data. */
static void build(struct osh_scoring_workspace **ws, struct osh_scoring_runtime *rt) {
    enum osh_status rc;

    write_detect_file();
    *ws = NULL;
    memset(rt, 0, sizeof(*rt));
    rc = osh_scoring_setup_from_path(DETECT_PATH, NULL, ws);
    ASSERT_TRUE(rc == OSH_OK);
    rc = osh_scoring_compile(*ws, NULL, rt);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(rt->noutputs == 2u);
    ASSERT_TRUE(rt->npages == 2u);
    rt->pages[0].acc.data[0] = 1.0;
    rt->pages[1].acc.data[0] = 2.0;
}

static void teardown(struct osh_scoring_workspace *ws, struct osh_scoring_runtime *rt) {
    osh_scoring_runtime_free(rt);
    osh_scoring_workspace_free(ws);
    remove(DETECT_PATH);
    remove(OUT0);
    remove(OUT1);
}

/* ---- osh_scoring_save_outputs validation branches ------------------------- */

static void test_save_outputs_errors(void) {
    struct osh_scoring_workspace *ws;
    struct osh_scoring_runtime rt;
    size_t const bad_idx[1] = {99u};
    size_t saved_noutputs;

    build(&ws, &rt);

    ASSERT_TRUE(osh_scoring_save_outputs(NULL, &rt, 5u, NULL, 0u) == OSH_EINVAL);
    ASSERT_TRUE(osh_scoring_save_outputs(ws, &rt, 0u, NULL, 0u) == OSH_EINVAL);    /* nstat must be > 0 */
    ASSERT_TRUE(osh_scoring_save_outputs(ws, &rt, 5u, bad_idx, 1u) == OSH_EINVAL); /* index out of range */

    /* Runtime with fewer outputs than the cold workspace -> ESTATE.  (A runtime
     * with *more* outputs is legitimate: multi-format blocks fan out, issue #308.) */
    saved_noutputs = rt.noutputs;
    rt.noutputs = saved_noutputs - 1u;
    ASSERT_TRUE(osh_scoring_save_outputs(ws, &rt, 5u, NULL, 0u) == OSH_ESTATE);
    rt.noutputs = saved_noutputs;

    teardown(ws, &rt);
}

/* ---- File sink + output selector through the snapshot path ----------------- */

static void test_file_sink_selector(void) {
    struct osh_scoring_workspace *ws;
    struct osh_scoring_runtime rt;
    struct osh_scoring_shadow shadow;
    struct osh_scoring_file_sink fs;
    struct osh_scoring_sink sink;
    size_t const want0[1] = {0u};

    build(&ws, &rt);
    remove(OUT0);
    remove(OUT1);

    ASSERT_TRUE(osh_scoring_file_sink_init(&fs, ws, &sink) == OSH_OK);
    ASSERT_TRUE(osh_scoring_shadow_init(&shadow, &rt) == OSH_OK);

    /* Selector: only output 0 is written. */
    ASSERT_TRUE(osh_scoring_snapshot_save(&sink, &shadow, 5u, want0, 1u) == OSH_OK);
    ASSERT_TRUE(file_exists(OUT0));
    ASSERT_TRUE(!file_exists(OUT1));

    /* NULL selector: all outputs are written. */
    ASSERT_TRUE(osh_scoring_snapshot_save(&sink, &shadow, 5u, NULL, 0u) == OSH_OK);
    ASSERT_TRUE(file_exists(OUT0));
    ASSERT_TRUE(file_exists(OUT1));

    /* file_sink_save rejects a missing workspace. */
    {
        struct osh_scoring_file_sink empty = {0};
        struct osh_scoring_sink bad_sink;
        bad_sink.save = sink.save;
        bad_sink.ctx = &empty; /* ws == NULL */
        ASSERT_TRUE(osh_scoring_snapshot_save(&bad_sink, &shadow, 5u, NULL, 0u) == OSH_EINVAL);
    }

    osh_scoring_shadow_free(&shadow);
    teardown(ws, &rt);
}

int main(void) {
    test_save_outputs_errors();
    test_file_sink_selector();
    printf("All osh_scoring_sink tests passed.\n");
    return 0;
}
