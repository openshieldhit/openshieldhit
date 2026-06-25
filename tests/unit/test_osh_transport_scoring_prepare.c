#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_app_osh.h"
#include "openshieldhit/scoring.h"
#include "openshieldhit/status.h"
#include "scoring/runtime/osh_scoring_compile.h"

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

    snprintf(path, path_cap, "osh_scoring_prepare_test_%d.tmp", tmp_counter++);
    fp = fopen(path, "w");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fputs(content, fp) >= 0);
    ASSERT_TRUE(fclose(fp) == 0);
}

static void test_compile_fixture_test01_detect(void) {
    char path[512];
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    enum osh_status rc;
    size_t page_idx;

    snprintf(path, sizeof(path), "%s/tests/fixtures/test01/detect.dat", OSH_PROJECT_SOURCE_DIR);

    rc = osh_scoring_setup_from_path(path, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(ws != NULL);

    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);

    ASSERT_TRUE(rt.nfilters == 2u);
    ASSERT_TRUE(rt.nsettings == 0u);
    ASSERT_TRUE(rt.ngeometries == 1u);
    ASSERT_TRUE(rt.noutputs == 2u);
    ASSERT_TRUE(rt.npages == 4u);

    ASSERT_TRUE(strcmp(rt.geometries[0].name, "MyMesh") == 0);
    ASSERT_TRUE(strcmp(rt.geometries[0].kind, "mesh") == 0);
    ASSERT_TRUE(rt.geometries[0].geo_kind == OSH_SCORING_GEO_MESH);
    ASSERT_TRUE(rt.geometries[0].nbins == 10u);
    ASSERT_TRUE(rt.geometries[0].first_page == 0u);
    ASSERT_TRUE(rt.geometries[0].npages == 4u);
    ASSERT_TRUE(rt.geometries[0].ngroups == 2u);
    ASSERT_TRUE(rt.geometries[0].groups[0].score_kind == OSH_SCORING_SCORE_ENERGY);
    ASSERT_TRUE(rt.geometries[0].groups[0].first_page == 0u);
    ASSERT_TRUE(rt.geometries[0].groups[0].npages == 2u);
    ASSERT_TRUE(rt.geometries[0].groups[1].score_kind == OSH_SCORING_SCORE_FLUENCE);
    ASSERT_TRUE(rt.geometries[0].groups[1].first_page == 2u);
    ASSERT_TRUE(rt.geometries[0].groups[1].npages == 2u);

    ASSERT_TRUE(strcmp(rt.outputs[0].filename, "NB_msh_energy.bdo") == 0);
    ASSERT_TRUE(strcmp(rt.outputs[0].fileformat, "BDO") == 0);
    ASSERT_TRUE(rt.outputs[0].geometry_idx == 0u);
    ASSERT_TRUE(rt.outputs[0].npages == 1u);

    page_idx = rt.outputs[0].page_indices[0];
    ASSERT_TRUE(page_idx < rt.npages);
    ASSERT_TRUE(strcmp(rt.pages[page_idx].quantity, "energy") == 0);
    ASSERT_TRUE(rt.pages[page_idx].score_kind == OSH_SCORING_SCORE_ENERGY);
    ASSERT_TRUE(rt.pages[page_idx].geometry_idx == 0u);
    ASSERT_TRUE(rt.pages[page_idx].output_idx == 0u);
    ASSERT_TRUE(rt.pages[page_idx].len == 10u);
    ASSERT_TRUE(rt.pages[page_idx].postproc == OSH_SCORING_POSTPROC_NORM);
    ASSERT_TRUE(rt.pages[page_idx].has_data2 == 0);
    ASSERT_TRUE(rt.pages[page_idx].divide == 0);
    ASSERT_TRUE(rt.pages[page_idx].acc.data != NULL);
    ASSERT_TRUE(rt.pages[page_idx].nflat_rules == 0u);

    ASSERT_TRUE(strcmp(rt.outputs[1].filename, "NB_msh_fluence.bdo") == 0);
    ASSERT_TRUE(strcmp(rt.outputs[1].fileformat, "BDO") == 0);
    ASSERT_TRUE(rt.outputs[1].geometry_idx == 0u);
    ASSERT_TRUE(rt.outputs[1].npages == 3u);
    ASSERT_TRUE(rt.outputs[1].page_indices[0] == 1u);
    ASSERT_TRUE(rt.outputs[1].page_indices[1] == 2u);
    ASSERT_TRUE(rt.outputs[1].page_indices[2] == 3u);

    page_idx = rt.outputs[1].page_indices[2];
    ASSERT_TRUE(page_idx < rt.npages);
    ASSERT_TRUE(strcmp(rt.pages[page_idx].quantity, "fluence") == 0);
    ASSERT_TRUE(rt.pages[page_idx].score_kind == OSH_SCORING_SCORE_FLUENCE);
    ASSERT_TRUE(rt.pages[page_idx].output_idx == 1u);
    ASSERT_TRUE(rt.pages[page_idx].geometry_idx == 0u);
    ASSERT_TRUE(rt.pages[page_idx].len == 10u);
    ASSERT_TRUE(rt.pages[page_idx].nflat_rules == 4u); /* MyFilter: 3 rules + Gen2: 1 rule */
    ASSERT_TRUE(rt.pages[page_idx].postproc == OSH_SCORING_POSTPROC_NORM);
    ASSERT_TRUE(rt.pages[page_idx].has_data2 == 0);
    ASSERT_TRUE(rt.pages[page_idx].divide == 0);
    ASSERT_TRUE(rt.pages[page_idx].acc.data != NULL);
    ASSERT_TRUE(rt.pages[page_idx].acc.data_var == NULL);
    ASSERT_TRUE(rt.pages[page_idx].acc.data2 == NULL);
    ASSERT_TRUE(rt.pages[page_idx].acc.data2_var == NULL);

    ASSERT_TRUE(rt.filters[0].nrules == 3u);
    ASSERT_TRUE(rt.filters[0].rules[0].field == OSH_SCORING_FILTER_FIELD_Z);
    ASSERT_TRUE(rt.filters[0].rules[0].op == OSH_SCORING_FILTER_OP_EQ);
    ASSERT_TRUE(rt.filters[0].rules[1].field == OSH_SCORING_FILTER_FIELD_A);
    ASSERT_TRUE(rt.filters[0].rules[1].op == OSH_SCORING_FILTER_OP_EQ);
    ASSERT_TRUE(rt.filters[0].rules[2].field == OSH_SCORING_FILTER_FIELD_E);
    ASSERT_TRUE(rt.filters[0].rules[2].op == OSH_SCORING_FILTER_OP_GT);

    ASSERT_TRUE(rt.filters[1].nrules == 1u);
    ASSERT_TRUE(rt.filters[1].rules[0].field == OSH_SCORING_FILTER_FIELD_GEN);
    ASSERT_TRUE(rt.filters[1].rules[0].op == OSH_SCORING_FILTER_OP_EQ);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
}

static void test_prepare_compiles_dose_page(void) {
    char path[512];
    char const *text = "Geometry Mesh\n"
                       "    Name MyMesh\n"
                       "    X 0.0 1.0 1\n"
                       "    Y 0.0 1.0 1\n"
                       "    Z 0.0 1.0 1\n"
                       "\n"
                       "Output\n"
                       "    Filename out.bdo\n"
                       "    Geo MyMesh\n"
                       "    Quantity DOSE\n";
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    enum osh_status rc;

    write_temp_file(path, sizeof(path), text);

    rc = osh_scoring_setup_from_path(path, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(ws != NULL);

    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(rt.npages == 1u);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    ASSERT_TRUE(remove(path) == 0);
}

int main(void) {
    test_compile_fixture_test01_detect();
    test_prepare_compiles_dose_page();
    return 0;
}
