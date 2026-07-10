/*
 * Unit tests for attaching several output formats to one scoring Output block
 * (issue #308).
 *
 * The headline guarantee: asking for more output formats must NOT grow the
 * scored-array memory.  One Output block scores a page-set once; each requested
 * format is a lightweight runtime output sharing the same pages.  These tests
 * lock that down at the level the guarantee lives — the memory estimate and the
 * compiled runtime — plus the surrounding filename/derivation and rejection
 * rules.
 *
 * Coverage:
 *   test_memory_matches_single    — estimate for TEXT+BDO == estimate for TEXT
 *   test_fanout_shares_pages      — one block, two formats: 2 outputs, 2 pages,
 *                                   both outputs point at the same pages
 *   test_derived_filenames        — stem + canonical extension per format
 *   test_single_format_verbatim   — one format keeps the filename byte-for-byte
 *   test_writes_all_files         — save emits one file per format
 *   test_collision_rejected       — two formats resolving to one path is an error
 *   test_rtdose_mixed_rejected    — RTDOSE cannot share a multi-format block
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_app_osh.h"
#include "openshieldhit/scoring.h"
#include "openshieldhit/status.h"
#include "scoring/runtime/osh_scoring_compile.h"
#include "scoring/save/osh_scoring_save.h"
#include "test_assert.h"

static int tmp_counter = 0;

static void write_temp_detect_with_filename(char *path,
                                            size_t path_cap,
                                            char const *filename,
                                            char const *fileformat_line) {
    FILE *fp;

    snprintf(path, path_cap, "osh_multiformat_test_%d.tmp", tmp_counter++);
    fp = fopen(path, "w");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fprintf(fp,
                        "Geometry Mesh\n"
                        "    Name MyMesh\n"
                        "    X 0.0 1.0 2\n"
                        "    Y 0.0 1.0 3\n"
                        "    Z 0.0 1.0 4\n"
                        "\n"
                        "Output\n"
                        "    Filename %s\n"
                        "    %s\n"
                        "    Geo MyMesh\n"
                        "    Quantity Energy\n"
                        "    Quantity Fluence\n",
                        filename,
                        fileformat_line)
                >= 0);
    ASSERT_TRUE(fclose(fp) == 0);
}

static void write_temp_detect(char *path, size_t path_cap, char const *fileformat_line) {
    write_temp_detect_with_filename(path, path_cap, "NB_msh", fileformat_line);
}

static int file_exists(char const *path) {
    FILE *fp = fopen(path, "r");
    if (fp) {
        (void) fclose(fp);
        return 1;
    }
    return 0;
}

/* The memory estimate reads the cold workspace, which holds one page-set per
 * Output block regardless of how many formats are attached — so TEXT+BDO must
 * estimate byte-for-byte the same accumulator memory as TEXT alone. */
static void test_memory_matches_single(void) {
    char path_multi[512];
    char path_single[512];
    struct osh_scoring_workspace *ws_multi = NULL;
    struct osh_scoring_workspace *ws_single = NULL;
    struct osh_scoring_mem_estimate est_multi;
    struct osh_scoring_mem_estimate est_single;

    write_temp_detect(path_multi, sizeof(path_multi), "FileFormat TEXT BDO");
    write_temp_detect(path_single, sizeof(path_single), "FileFormat TEXT");

    ASSERT_TRUE(osh_scoring_setup_from_path(path_multi, NULL, &ws_multi) == OSH_OK);
    ASSERT_TRUE(osh_scoring_setup_from_path(path_single, NULL, &ws_single) == OSH_OK);

    ASSERT_TRUE(osh_scoring_estimate_memory(ws_multi, &est_multi) == OSH_OK);
    ASSERT_TRUE(osh_scoring_estimate_memory(ws_single, &est_single) == OSH_OK);

    ASSERT_TRUE(est_multi.accum_bytes == est_single.accum_bytes);
    ASSERT_TRUE(est_multi.npages == est_single.npages);
    ASSERT_TRUE(est_multi.shadow_bytes == est_single.shadow_bytes);
    ASSERT_TRUE(est_multi.largest_page_bytes == est_single.largest_page_bytes);
    /* 2 quantities × 24 bins × 8 bytes, counted once. */
    ASSERT_TRUE(est_multi.npages == 2u);
    ASSERT_TRUE(est_multi.accum_bytes == (size_t) 2u * 24u * sizeof(double));

    osh_scoring_workspace_free(ws_multi);
    osh_scoring_workspace_free(ws_single);
    remove(path_multi);
    remove(path_single);
}

/* One block with two formats compiles into two runtime outputs that share the
 * same pages: 2 pages total (not 4), and identical page-index lists in separate
 * arrays — proving the accumulators are shared, not duplicated. */
static void test_fanout_shares_pages(void) {
    char path[512];
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    size_t j;

    write_temp_detect(path, sizeof(path), "FileFormat TEXT BDO");
    ASSERT_TRUE(osh_scoring_setup_from_path(path, NULL, &ws) == OSH_OK);

    memset(&rt, 0, sizeof(rt));
    ASSERT_TRUE(osh_scoring_compile(ws, NULL, &rt) == OSH_OK);

    /* Two formats -> two runtime outputs, but still only two scored pages. */
    ASSERT_TRUE(rt.noutputs == 2u);
    ASSERT_TRUE(rt.npages == 2u);

    /* Primary keeps the first-listed format; the extra carries the second. */
    ASSERT_TRUE(strcmp(rt.outputs[0].fileformat, "text") == 0);
    ASSERT_TRUE(strcmp(rt.outputs[1].fileformat, "bdo") == 0);
    ASSERT_TRUE(rt.outputs[0].npages == 2u);
    ASSERT_TRUE(rt.outputs[1].npages == 2u);

    /* Same pages, different arrays: shared indices, not shared storage handle. */
    ASSERT_TRUE(rt.outputs[0].page_indices != rt.outputs[1].page_indices);
    for (j = 0u; j < rt.outputs[0].npages; ++j) {
        ASSERT_TRUE(rt.outputs[0].page_indices[j] == rt.outputs[1].page_indices[j]);
        ASSERT_TRUE(rt.outputs[0].page_indices[j] < rt.npages);
    }
    ASSERT_TRUE(rt.outputs[0].geometry_idx == rt.outputs[1].geometry_idx);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(path);
}

/* With more than one format the Filename is a stem; each format gets its
 * canonical extension. */
static void test_derived_filenames(void) {
    char path[512];
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;

    write_temp_detect_with_filename(path, sizeof(path), "NB_msh.DAT", "FileFormat TEXT BDO");
    ASSERT_TRUE(osh_scoring_setup_from_path(path, NULL, &ws) == OSH_OK);

    memset(&rt, 0, sizeof(rt));
    ASSERT_TRUE(osh_scoring_compile(ws, NULL, &rt) == OSH_OK);

    ASSERT_TRUE(strcmp(rt.outputs[0].filename, "NB_msh.dat") == 0);
    ASSERT_TRUE(strcmp(rt.outputs[1].filename, "NB_msh.bdo") == 0);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(path);
}

/* A single format leaves the filename exactly as written — no extension is
 * appended or stripped — so every existing detect.dat is unaffected. */
static void test_single_format_verbatim(void) {
    char path[512];
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;

    write_temp_detect(path, sizeof(path), "FileFormat TEXT");
    ASSERT_TRUE(osh_scoring_setup_from_path(path, NULL, &ws) == OSH_OK);

    memset(&rt, 0, sizeof(rt));
    ASSERT_TRUE(osh_scoring_compile(ws, NULL, &rt) == OSH_OK);

    ASSERT_TRUE(rt.noutputs == 1u);
    ASSERT_TRUE(strcmp(rt.outputs[0].filename, "NB_msh") == 0);
    ASSERT_TRUE(strcmp(rt.outputs[0].fileformat, "text") == 0);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(path);
}

/* Saving a two-format block emits one file per format. */
static void test_writes_all_files(void) {
    char path[512];
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;

    write_temp_detect(path, sizeof(path), "FileFormat TEXT BDO");
    ASSERT_TRUE(osh_scoring_setup_from_path(path, NULL, &ws) == OSH_OK);

    memset(&rt, 0, sizeof(rt));
    ASSERT_TRUE(osh_scoring_compile(ws, NULL, &rt) == OSH_OK);

    rt.pages[0].acc.data[0] = 1.0;
    rt.pages[1].acc.data[0] = 2.0;

    remove("NB_msh.dat");
    remove("NB_msh.bdo");
    ASSERT_TRUE(osh_scoring_save(ws, &rt, 10u) == OSH_OK);
    ASSERT_TRUE(file_exists("NB_msh.dat"));
    ASSERT_TRUE(file_exists("NB_msh.bdo"));

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(path);
    remove("NB_msh.dat");
    remove("NB_msh.bdo");
}

/* Two formats that canonicalise to the same path (TEXT and DAT both -> .dat)
 * are rejected at compile time. */
static void test_collision_rejected(void) {
    char path[512];
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;

    write_temp_detect(path, sizeof(path), "FileFormat TEXT DAT");
    ASSERT_TRUE(osh_scoring_setup_from_path(path, NULL, &ws) == OSH_OK);

    memset(&rt, 0, sizeof(rt));
    ASSERT_TRUE(osh_scoring_compile(ws, NULL, &rt) == OSH_EINVAL);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(path);
}

/* RTDOSE needs exactly one dose page and cannot share a multi-page block, so
 * combining it with another format in one block is rejected. */
static void test_rtdose_mixed_rejected(void) {
    char path[512];
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;

    write_temp_detect(path, sizeof(path), "FileFormat BDO RTDOSE");
    ASSERT_TRUE(osh_scoring_setup_from_path(path, NULL, &ws) == OSH_OK);

    memset(&rt, 0, sizeof(rt));
    ASSERT_TRUE(osh_scoring_compile(ws, NULL, &rt) == OSH_ENOTSUP);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(path);
}

int main(void) {
    test_memory_matches_single();
    test_fanout_shares_pages();
    test_derived_filenames();
    test_single_format_verbatim();
    test_writes_all_files();
    test_collision_rejected();
    test_rtdose_mixed_rejected();
    return 0;
}
