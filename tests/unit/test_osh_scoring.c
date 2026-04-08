#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
    rc = osh_scoring_setup_from_path(path, &ws);
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
    ASSERT_TRUE(strcmp(geo->kind, "Mesh") == 0);

    out0 = osh_scoring_output_by_filename(ws, "NB_msh_energy.bdo");
    ASSERT_TRUE(out0 != NULL);
    ASSERT_TRUE(strcmp(out0->geometry_name, "MyMesh") == 0);
    ASSERT_TRUE(out0->npages == 1u);
    ASSERT_TRUE(strcmp(out0->pages[0].quantity, "ENERGY") == 0);
    ASSERT_TRUE(out0->pages[0].nfilter_names == 0u);

    out1 = osh_scoring_output_by_filename(ws, "NB_msh_fluence.bdo");
    ASSERT_TRUE(out1 != NULL);
    ASSERT_TRUE(strcmp(out1->geometry_name, "MyMesh") == 0);
    ASSERT_TRUE(out1->npages == 3u);
    ASSERT_TRUE(strcmp(out1->pages[0].quantity, "ENERGY") == 0);
    ASSERT_TRUE(strcmp(out1->pages[1].quantity, "FLUENCE") == 0);
    ASSERT_TRUE(strcmp(out1->pages[2].quantity, "FLUENCE") == 0);
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
    rc = osh_scoring_setup_from_path(path, &ws);
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
    rc = osh_scoring_setup_from_path(path, &ws);
    ASSERT_TRUE(rc == OSH_EPARSE);
    ASSERT_TRUE(ws == NULL);
    remove(path);
}

int main(void) {
    test_parse_fixture_test01_detect();
    test_parse_settings_section();
    test_reject_output_missing_filename();
    return 0;
}
