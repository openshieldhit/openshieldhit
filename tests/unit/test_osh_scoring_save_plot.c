/*
 * Native SVG plot writer (issue #238 prototype).
 *
 * The feature is compiled in only under -DOSH_ENABLE_PLOT=ON, which also defines
 * the OSH_ENABLE_PLOT macro globally.  This test drives only the public
 * osh_scoring_save() dispatcher (never the gated writer symbol directly), so it
 * links and runs in both configurations:
 *
 *   - built with the feature: a supported shape yields OSH_OK and a real SVG doc;
 *   - stock build (default):  `FileFormat SVG` is unrecognised -> OSH_ENOTSUP,
 *     proving the default binary stays plotting-free.
 *
 * Every case runs osh_scoring_postprocess() before saving, matching the real
 * save pipeline (and the writer's contract, which expects postprocessed
 * accumulators — differential bin-width division and volume normalisation both
 * happen there).
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
#ifdef OSH_ENABLE_PLOT
static int file_contains(char const *path, char const *needle);
#endif

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

#ifdef OSH_ENABLE_PLOT
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(file_contains("out_plot_profile.svg", "<svg"));
    ASSERT_TRUE(file_contains("out_plot_profile.svg", "<polyline"));
    ASSERT_TRUE(file_contains("out_plot_profile.svg", "ENERGY"));
    ASSERT_TRUE(file_contains("out_plot_profile.svg", "Z [cm]"));
#else
    /* Stock build: the SVG format is not compiled in, so it is unsupported. */
    ASSERT_TRUE(rc == OSH_ENOTSUP);
#endif

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

#ifdef OSH_ENABLE_PLOT
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(file_contains("out_plot_spectrum.svg", "<polyline"));
    ASSERT_TRUE(file_contains("out_plot_spectrum.svg", "E [MeV]"));   /* differential x-axis */
    ASSERT_TRUE(file_contains("out_plot_spectrum.svg", "/cm^2/MeV")); /* differential y unit */
#else
    ASSERT_TRUE(rc == OSH_ENOTSUP);
#endif

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

#ifdef OSH_ENABLE_PLOT
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(file_contains("out_plot_zspectrum.svg", "<polyline"));
    ASSERT_TRUE(file_contains("out_plot_zspectrum.svg", "E [MeV]"));
#else
    ASSERT_TRUE(rc == OSH_ENOTSUP);
#endif

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(DETECT_PATH);
    remove("out_plot_zspectrum.svg");
}

/* A 2-D mesh (two non-singleton axes) is out of scope for the 1-D prototype and
 * must be declined cleanly.  It is OSH_ENOTSUP either way: the writer rejects the
 * shape when built in, the dispatcher rejects the format when not. */
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

#ifdef OSH_ENABLE_PLOT
static int file_contains(char const *path, char const *needle) {
    FILE *fp;
    char line[1024];
    int found;

    fp = fopen(path, "r");
    if (!fp) {
        return 0;
    }
    found = 0;
    while (fgets(line, sizeof(line), fp) != NULL) {
        if (strstr(line, needle) != NULL) {
            found = 1;
            break;
        }
    }
    fclose(fp);
    return found;
}
#endif
