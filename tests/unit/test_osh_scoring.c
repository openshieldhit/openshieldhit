#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_app_osh.h"
#include "common/osh_rc.h"
#include "scoring/osh_scoring.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

static int tmp_counter = 0;

static void write_temp_file(char *path, size_t path_cap, char const *content) {
    FILE *fp;

    snprintf(path, path_cap, "osh_scoring_test_%d.tmp", tmp_counter++);
    fp = fopen(path, "w");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fputs(content, fp) >= 0);
    ASSERT_TRUE(fclose(fp) == 0);
}

static void test_parse_fixture_test01_detect(void) {
    char path[512];
    struct osh_scoring_workspace *ws;
    struct osh_scoring_geometry_def const *geo;
    struct osh_scoring_output_def const *out0;
    struct osh_scoring_output_def const *out1;
    enum osh_status rc;

    snprintf(path, sizeof(path), "%s/tests/fixtures/test01/detect.dat", OSH_PROJECT_SOURCE_DIR);

    ws = NULL;
    rc = osh_scoring_setup_from_path(path, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(ws != NULL);

    ASSERT_TRUE(ws->nfilters == 2u);
    ASSERT_TRUE(ws->nsettings == 0u);
    ASSERT_TRUE(ws->ngeometries == 1u);
    ASSERT_TRUE(ws->noutputs == 2u);

    ASSERT_TRUE(osh_scoring_filter_by_name(ws, "MyFilter") != NULL);
    ASSERT_TRUE(osh_scoring_filter_by_name(ws, "Gen2") != NULL);

    geo = osh_scoring_geometry_by_name(ws, "MyMesh");
    ASSERT_TRUE(geo != NULL);
    ASSERT_TRUE(strcmp(geo->kind, "mesh") == 0);

    out0 = osh_scoring_output_by_filename(ws, "NB_msh_energy.bdo");
    ASSERT_TRUE(out0 != NULL);
    ASSERT_TRUE(strcmp(out0->geometry_name, "MyMesh") == 0);
    ASSERT_TRUE(out0->npages == 1u);
    ASSERT_TRUE(strcmp(out0->pages[0].quantity, "energy") == 0);
    ASSERT_TRUE(out0->pages[0].nfilter_names == 0u);

    out1 = osh_scoring_output_by_filename(ws, "NB_msh_fluence.bdo");
    ASSERT_TRUE(out1 != NULL);
    ASSERT_TRUE(strcmp(out1->geometry_name, "MyMesh") == 0);
    ASSERT_TRUE(out1->npages == 3u);
    ASSERT_TRUE(strcmp(out1->pages[0].quantity, "energy") == 0);
    ASSERT_TRUE(strcmp(out1->pages[1].quantity, "fluence") == 0);
    ASSERT_TRUE(strcmp(out1->pages[2].quantity, "fluence") == 0);
    ASSERT_TRUE(out1->pages[2].nfilter_names == 2u);
    ASSERT_TRUE(strcmp(out1->pages[2].filter_names[0], "MyFilter") == 0);
    ASSERT_TRUE(strcmp(out1->pages[2].filter_names[1], "Gen2") == 0);

    osh_scoring_workspace_free(ws);
}

static void test_parse_settings_section(void) {
    char path[512];
    char const *text = "Settings\n"
                       "    Name DoseWater\n"
                       "    Rescale 2.5\n"
                       "    Offset 1.25\n"
                       "    Density 0.998\n"
                       "    SiteDiam 0.75\n"
                       "    MaxCount 123\n"
                       "    Medium 7\n"
                       "    NKMedium 3\n"
                       "\n"
                       "Geometry Mesh\n"
                       "    Name MyMesh\n"
                       "\n"
                       "Output\n"
                       "    Filename out.bdo\n"
                       "    Geo MyMesh\n"
                       "    Quantity ENERGY\n";
    struct osh_scoring_workspace *ws;
    struct osh_scoring_settings_def const *set;
    enum osh_status rc;

    write_temp_file(path, sizeof(path), text);

    ws = NULL;
    rc = osh_scoring_setup_from_path(path, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(ws != NULL);
    ASSERT_TRUE(ws->nsettings == 1u);

    set = osh_scoring_settings_by_name(ws, "DoseWater");
    ASSERT_TRUE(set != NULL);
    ASSERT_TRUE(set->has_rescale);
    ASSERT_TRUE(set->has_offset);
    ASSERT_TRUE(set->has_density_g_cm3);
    ASSERT_TRUE(set->has_site_diameter_um);
    ASSERT_TRUE(set->has_npart);
    ASSERT_TRUE(set->has_medium);
    ASSERT_TRUE(set->has_nkmedium);
    ASSERT_TRUE(set->rescale == 2.5);
    ASSERT_TRUE(set->offset == 1.25);
    ASSERT_TRUE(set->density_g_cm3 == 0.998);
    ASSERT_TRUE(set->site_diameter_um == 0.75);
    ASSERT_TRUE(set->npart == 123u);
    ASSERT_TRUE(set->medium == 7);
    ASSERT_TRUE(set->nkmedium == 3);

    osh_scoring_workspace_free(ws);
    remove(path);
}

static void test_reject_output_missing_filename(void) {
    char path[512];
    char const *text = "Geometry Mesh\n"
                       "    Name MyMesh\n"
                       "\n"
                       "Output\n"
                       "    Geo MyMesh\n"
                       "    Quantity ENERGY\n";
    struct osh_scoring_workspace *ws;
    enum osh_status rc;

    write_temp_file(path, sizeof(path), text);

    ws = NULL;
    rc = osh_scoring_setup_from_path(path, NULL, &ws);
    ASSERT_TRUE(rc == OSH_EPARSE);
    ASSERT_TRUE(ws == NULL);
    remove(path);
}

static void test_parse_filter_rules(void) {
    char path[512];
    char const *text = "Filter\n"
                       "    Name MyFilter\n"
                       "    Z = 6\n"
                       "    A = 12\n"
                       "    E > 0.1\n"
                       "\n"
                       "Geometry Mesh\n"
                       "    Name G\n"
                       "\n"
                       "Output\n"
                       "    Filename out.bdo\n"
                       "    Geo G\n"
                       "    Quantity DOSE\n";
    struct osh_scoring_workspace *ws;
    struct osh_scoring_filter_def const *fil;
    enum osh_status rc;

    write_temp_file(path, sizeof(path), text);
    ws = NULL;
    rc = osh_scoring_setup_from_path(path, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(ws != NULL);

    fil = osh_scoring_filter_by_name(ws, "MyFilter");
    ASSERT_TRUE(fil != NULL);
    ASSERT_TRUE(fil->nrules == 3u);
    ASSERT_TRUE(strcmp(fil->rules[0].field, "Z") == 0);
    ASSERT_TRUE(strcmp(fil->rules[0].op, "=") == 0);
    ASSERT_TRUE(fil->rules[0].value == 6.0);
    ASSERT_TRUE(strcmp(fil->rules[1].field, "A") == 0);
    ASSERT_TRUE(fil->rules[1].value == 12.0);
    ASSERT_TRUE(strcmp(fil->rules[2].field, "E") == 0);
    ASSERT_TRUE(strcmp(fil->rules[2].op, ">") == 0);
    ASSERT_TRUE(fil->rules[2].value == 0.1);

    osh_scoring_workspace_free(ws);
    remove(path);
}

static void test_parse_geometry_axes(void) {
    char path[512];
    char const *text = "Geometry Mesh\n"
                       "    Name MyMesh\n"
                       "    X -5.0  5.0  10\n"
                       "    Y -5.0  5.0  10\n"
                       "    Z  0.0 30.0  60\n"
                       "\n"
                       "Geometry Cyl\n"
                       "    Name MyCyl\n"
                       "    R  0.0  3.0  5\n"
                       "    Z  0.0 30.0 60\n"
                       "\n"
                       "Output\n"
                       "    Filename mesh.bdo\n"
                       "    Geo MyMesh\n"
                       "    Quantity DOSE\n"
                       "\n"
                       "Output\n"
                       "    Filename cyl.bdo\n"
                       "    Geo MyCyl\n"
                       "    Quantity DOSE\n";
    struct osh_scoring_workspace *ws;
    struct osh_scoring_geometry_def const *mesh;
    struct osh_scoring_geometry_def const *cyl;
    enum osh_status rc;

    write_temp_file(path, sizeof(path), text);
    ws = NULL;
    rc = osh_scoring_setup_from_path(path, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(ws != NULL);
    ASSERT_TRUE(ws->ngeometries == 2u);

    mesh = osh_scoring_geometry_by_name(ws, "MyMesh");
    ASSERT_TRUE(mesh != NULL);
    ASSERT_TRUE(strcmp(mesh->kind, "mesh") == 0);
    ASSERT_TRUE(mesh->naxes == 3u);
    ASSERT_TRUE(strcmp(mesh->axes[0].label, "X") == 0);
    ASSERT_TRUE(mesh->axes[0].lo == -5.0);
    ASSERT_TRUE(mesh->axes[0].hi == 5.0);
    ASSERT_TRUE(mesh->axes[0].nbins == 10);
    ASSERT_TRUE(strcmp(mesh->axes[2].label, "Z") == 0);
    ASSERT_TRUE(mesh->axes[2].nbins == 60);

    cyl = osh_scoring_geometry_by_name(ws, "MyCyl");
    ASSERT_TRUE(cyl != NULL);
    ASSERT_TRUE(strcmp(cyl->kind, "cyl") == 0);
    ASSERT_TRUE(cyl->naxes == 2u);
    ASSERT_TRUE(strcmp(cyl->axes[0].label, "R") == 0);
    ASSERT_TRUE(cyl->axes[0].hi == 3.0);
    ASSERT_TRUE(cyl->axes[0].nbins == 5);

    osh_scoring_workspace_free(ws);
    remove(path);
}

static void test_fixture_test01_filter_rules(void) {
    char path[512];
    struct osh_scoring_workspace *ws;
    struct osh_scoring_filter_def const *mf;
    struct osh_scoring_filter_def const *g2;
    enum osh_status rc;

    snprintf(path, sizeof(path), "%s/tests/fixtures/test01/detect.dat", OSH_PROJECT_SOURCE_DIR);
    ws = NULL;
    rc = osh_scoring_setup_from_path(path, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);

    mf = osh_scoring_filter_by_name(ws, "MyFilter");
    ASSERT_TRUE(mf != NULL);
    ASSERT_TRUE(mf->nrules == 3u); /* Z=6, A=12, E>0.1 */

    g2 = osh_scoring_filter_by_name(ws, "Gen2");
    ASSERT_TRUE(g2 != NULL);
    ASSERT_TRUE(g2->nrules == 1u); /* GEN=0 */
    ASSERT_TRUE(strcmp(g2->rules[0].field, "GEN") == 0);
    ASSERT_TRUE(g2->rules[0].value == 0.0);

    /* Geometry axes */
    {
        struct osh_scoring_geometry_def const *geo = osh_scoring_geometry_by_name(ws, "MyMesh");
        ASSERT_TRUE(geo != NULL);
        ASSERT_TRUE(geo->naxes == 3u);
        ASSERT_TRUE(strcmp(geo->axes[0].label, "X") == 0);
        ASSERT_TRUE(geo->axes[0].lo == -0.5);
        ASSERT_TRUE(geo->axes[0].hi == 0.5);
        ASSERT_TRUE(geo->axes[0].nbins == 1);
        ASSERT_TRUE(strcmp(geo->axes[2].label, "Z") == 0);
        ASSERT_TRUE(geo->axes[2].nbins == 10);
    }

    osh_scoring_workspace_free(ws);
}

int main(void) {
    test_parse_fixture_test01_detect();
    test_parse_settings_section();
    test_reject_output_missing_filename();
    test_parse_filter_rules();
    test_parse_geometry_axes();
    test_fixture_test01_filter_rules();
    return 0;
}
