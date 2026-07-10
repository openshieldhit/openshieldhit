/*
 * Native SVG plot writer (issue #238).
 *
 * SVG output is always compiled in, reached through the public
 * osh_scoring_save() dispatcher via `FileFormat SVG`.  Every case runs
 * osh_scoring_postprocess() before saving, matching the real save pipeline (the
 * writer's contract expects postprocessed accumulators — differential bin-width
 * division and volume normalisation both happen there).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_app_osh.h"
#include "openshieldhit/scoring.h"
#include "openshieldhit/status.h"
#include "scoring/runtime/osh_scoring_compile.h"
#include "scoring/runtime/osh_scoring_postprocess.h"
#include "scoring/save/osh_scoring_save.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

#define DETECT_PATH "osh_scoring_save_plot_detect.tmp"

static void write_detect_file(char const *content);
static int file_contains(char const *path, char const *needle);
static int file_count(char const *path, char const *needle);

/* A 1-D depth profile (1 x 1 x N mesh) is the flagship case: one polyline. */
static void test_plot_svg_mesh_profile(void) {
    char const *detect_text = "Geometry Mesh\n"
                              "    Name G\n"
                              "    X 0 1 1\n"
                              "    Y 0 1 1\n"
                              "    Z 0 8 8\n"
                              "\n"
                              "Output\n"
                              "    Filename out_plot_profile.svg\n"
                              "    FileFormat SVG\n"
                              "    Geo G\n"
                              "    Quantity ENERGY\n";
    struct osh_scoring_workspace *ws;
    struct osh_scoring_runtime rt;
    enum osh_status rc;
    size_t i;

    write_detect_file(detect_text);

    ws = NULL;
    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_setup_from_path(DETECT_PATH, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(rt.npages == 1u);

    /* A crude Bragg-like bump so the polyline is non-degenerate. */
    for (i = 0; i < rt.pages[0].len; ++i) {
        rt.pages[0].acc.data[i] = 1.0 + (double) (i == 5u ? 9 : 0);
    }
    rc = osh_scoring_postprocess(&rt);
    ASSERT_TRUE(rc == OSH_OK);

    rc = osh_scoring_save(ws, &rt, 10u);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(file_contains("out_plot_profile.svg", "<svg"));
    ASSERT_TRUE(file_contains("out_plot_profile.svg", "<polyline"));
    ASSERT_TRUE(file_contains("out_plot_profile.svg", "ENERGY"));
    ASSERT_TRUE(file_contains("out_plot_profile.svg", "Z [cm]"));

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(DETECT_PATH);
    remove("out_plot_profile.svg");
}

/* An energy spectrum is a differential page over a single voxel (1 x 1 x 1
 * mesh).  Also a 1-D object: x is the differential (energy) axis, log-scaled. */
static void test_plot_svg_mesh_spectrum(void) {
    char const *detect_text = "Geometry Mesh\n"
                              "    Name G\n"
                              "    X 0 1 1\n"
                              "    Y 0 1 1\n"
                              "    Z 0 1 1\n"
                              "\n"
                              "Output\n"
                              "    Filename out_plot_spectrum.svg\n"
                              "    FileFormat SVG\n"
                              "    Geo G\n"
                              "    Quantity Fluence\n"
                              "    Diff1 1 100 16 LOG\n"
                              "    Diff1Type EKIN\n";
    struct osh_scoring_workspace *ws;
    struct osh_scoring_runtime rt;
    enum osh_status rc;
    size_t i;

    write_detect_file(detect_text);

    ws = NULL;
    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_setup_from_path(DETECT_PATH, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(rt.npages == 1u);
    ASSERT_TRUE(rt.pages[0].diff_nbins == 16u);
    ASSERT_TRUE(rt.pages[0].diff_stride == 1u); /* single spatial bin */

    for (i = 0; i < rt.pages[0].len; ++i) {
        rt.pages[0].acc.data[i] = 1.0 + (double) i;
    }
    rc = osh_scoring_postprocess(&rt);
    ASSERT_TRUE(rc == OSH_OK);

    rc = osh_scoring_save(ws, &rt, 10u);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(file_contains("out_plot_spectrum.svg", "<polyline"));
    ASSERT_TRUE(file_contains("out_plot_spectrum.svg", "E [MeV]"));   /* differential x-axis */
    ASSERT_TRUE(file_contains("out_plot_spectrum.svg", "/cm^2/MeV")); /* differential y unit */

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(DETECT_PATH);
    remove("out_plot_spectrum.svg");
}

/* The same spectrum scored over a single Zone.  A one-zone geometry is also a
 * 0-D spatial bin (diff_stride == 1), so it takes the same spectrum path. */
static void test_plot_svg_zone_spectrum(void) {
    char const *detect_text = "Geometry Zone\n"
                              "    Name Z\n"
                              "    Zone Target\n"
                              "    Volume 1.0\n"
                              "\n"
                              "Output\n"
                              "    Filename out_plot_zspectrum.svg\n"
                              "    FileFormat SVG\n"
                              "    Geo Z\n"
                              "    Quantity Fluence\n"
                              "    Diff1 1 100 8 LOG\n"
                              "    Diff1Type EKIN\n";
    struct osh_scoring_workspace *ws;
    struct osh_scoring_runtime rt;
    enum osh_status rc;
    size_t i;

    write_detect_file(detect_text);

    ws = NULL;
    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_setup_from_path(DETECT_PATH, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(ws != NULL);
    ASSERT_TRUE(ws->geometries[0].nzone_indices == 1u);

    /* Resolve the single Zone selector to a transport zone id (as the ASCII/BDO
     * zone save test does) so compile can proceed without a geo.dat workspace. */
    ws->geometries[0].zone_indices = (size_t *) calloc(1u, sizeof(*ws->geometries[0].zone_indices));
    ASSERT_TRUE(ws->geometries[0].zone_indices != NULL);
    ws->geometries[0].zone_indices[0] = 3u;

    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(rt.npages == 1u);
    ASSERT_TRUE(rt.pages[0].diff_stride == 1u); /* single zone */

    for (i = 0; i < rt.pages[0].len; ++i) {
        rt.pages[0].acc.data[i] = 1.0 + (double) i;
    }
    rc = osh_scoring_postprocess(&rt);
    ASSERT_TRUE(rc == OSH_OK);

    rc = osh_scoring_save(ws, &rt, 10u);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(file_contains("out_plot_zspectrum.svg", "<polyline"));
    ASSERT_TRUE(file_contains("out_plot_zspectrum.svg", "E [MeV]"));

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(DETECT_PATH);
    remove("out_plot_zspectrum.svg");
}

/* Mixed output on a 1-D depth mesh: two plain (1-D profile) pages plus one page
 * that adds a Diff1 axis (making it 2-D spatial x energy on this geometry).
 * The SVG must plot only the two truly-1-D pages and skip the differential one:
 * exactly two polylines, a spatial x-axis, and no differential x-axis. */
static void test_plot_svg_mixed_pages_profile(void) {
    char const *detect_text = "Geometry Mesh\n"
                              "    Name G\n"
                              "    X 0 1 1\n"
                              "    Y 0 1 1\n"
                              "    Z 0 8 8\n"
                              "\n"
                              "Output\n"
                              "    Filename out_plot_mixed.svg\n"
                              "    FileFormat SVG\n"
                              "    Geo G\n"
                              "    Quantity Dose\n"
                              "    Quantity Fluence\n"
                              "    Quantity Fluence\n"
                              "    Diff1 1 100 8 LOG\n"
                              "    Diff1Type EKIN\n";
    struct osh_scoring_workspace *ws;
    struct osh_scoring_runtime rt;
    enum osh_status rc;
    size_t p;
    size_t i;

    write_detect_file(detect_text);

    ws = NULL;
    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_setup_from_path(DETECT_PATH, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(rt.npages == 3u); /* Dose, Fluence, Fluence+Diff1 */

    for (p = 0; p < rt.npages; ++p) {
        for (i = 0; i < rt.pages[p].len; ++i) {
            rt.pages[p].acc.data[i] = 1.0 + (double) i;
        }
    }
    rc = osh_scoring_postprocess(&rt);
    ASSERT_TRUE(rc == OSH_OK);

    rc = osh_scoring_save(ws, &rt, 10u);
    ASSERT_TRUE(rc == OSH_OK);
    /* Only the two plain pages are 1-D here; the differential page is skipped. */
    ASSERT_TRUE(file_count("out_plot_mixed.svg", "<polyline") == 2);
    ASSERT_TRUE(file_contains("out_plot_mixed.svg", "Z [cm]"));   /* spatial profile */
    ASSERT_TRUE(!file_contains("out_plot_mixed.svg", "E [MeV]")); /* not a spectrum */
    ASSERT_TRUE(file_contains("out_plot_mixed.svg", "DOSE"));
    ASSERT_TRUE(file_contains("out_plot_mixed.svg", "FLUENCE"));

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(DETECT_PATH);
    remove("out_plot_mixed.svg");
}

/* A 2-D mesh (two non-singleton axes) fits neither shape and is declined. */
static void test_plot_svg_rejects_2d(void) {
    char const *detect_text = "Geometry Mesh\n"
                              "    Name G\n"
                              "    X 0 2 2\n"
                              "    Y 0 1 1\n"
                              "    Z 0 2 2\n"
                              "\n"
                              "Output\n"
                              "    Filename out_plot_2d.svg\n"
                              "    FileFormat SVG\n"
                              "    Geo G\n"
                              "    Quantity ENERGY\n";
    struct osh_scoring_workspace *ws;
    struct osh_scoring_runtime rt;
    enum osh_status rc;

    write_detect_file(detect_text);

    ws = NULL;
    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_setup_from_path(DETECT_PATH, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);

    rt.pages[0].acc.data[0] = 1.0;
    rc = osh_scoring_postprocess(&rt);
    ASSERT_TRUE(rc == OSH_OK);
    rc = osh_scoring_save(ws, &rt, 1u);
    ASSERT_TRUE(rc == OSH_ENOTSUP);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(DETECT_PATH);
    remove("out_plot_2d.svg");
}

int main(void) {
    test_plot_svg_mesh_profile();
    test_plot_svg_mesh_spectrum();
    test_plot_svg_zone_spectrum();
    test_plot_svg_mixed_pages_profile();
    test_plot_svg_rejects_2d();
    return 0;
}

static void write_detect_file(char const *content) {
    FILE *fp;

    fp = fopen(DETECT_PATH, "w");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fputs(content, fp) >= 0);
    ASSERT_TRUE(fclose(fp) == 0);
}

static int file_contains(char const *path, char const *needle) {
    return file_count(path, needle) > 0;
}

/* Count occurrences of @p needle across the whole file. */
static int file_count(char const *path, char const *needle) {
    FILE *fp;
    char line[4096];
    int n;

    fp = fopen(path, "r");
    if (!fp) {
        return -1;
    }
    n = 0;
    while (fgets(line, sizeof(line), fp) != NULL) {
        char const *pos = line;
        while ((pos = strstr(pos, needle)) != NULL) {
            ++n;
            pos += 1;
        }
    }
    fclose(fp);
    return n;
}
