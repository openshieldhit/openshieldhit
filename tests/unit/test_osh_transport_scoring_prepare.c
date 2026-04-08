#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_rc.h"
#include "scoring/osh_scoring.h"
#include "transport/prepare/osh_transport_scoring_prepare.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

static void test_compile_fixture_test01_detect(void) {
    char path[512];
    struct osh_scoring_workspace *ws = NULL;
    struct osh_transport_scoring_runtime rt;
    enum osh_status rc;
    size_t page_idx;

    snprintf(path, sizeof(path), "%s/tests/fixtures/test01/detect.dat", OSH_PROJECT_SOURCE_DIR);

    rc = osh_scoring_setup_from_path(path, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(ws != NULL);

    memset(&rt, 0, sizeof(rt));
    rc = osh_transport_scoring_prepare(ws, &rt);
    ASSERT_TRUE(rc == OSH_OK);

    ASSERT_TRUE(rt.nfilters == 2u);
    ASSERT_TRUE(rt.nsettings == 0u);
    ASSERT_TRUE(rt.ngeometries == 1u);
    ASSERT_TRUE(rt.noutputs == 2u);
    ASSERT_TRUE(rt.npages == 4u);

    ASSERT_TRUE(strcmp(rt.geometries[0].name, "MyMesh") == 0);
    ASSERT_TRUE(strcmp(rt.geometries[0].kind, "Mesh") == 0);
    ASSERT_TRUE(rt.geometries[0].nbins == 10u);
    ASSERT_TRUE(rt.geometries[0].first_page == 0u);
    ASSERT_TRUE(rt.geometries[0].npages == 4u);

    ASSERT_TRUE(strcmp(rt.outputs[0].filename, "NB_msh_energy.bdo") == 0);
    ASSERT_TRUE(strcmp(rt.outputs[0].fileformat, "BDO") == 0);
    ASSERT_TRUE(rt.outputs[0].geometry_idx == 0u);
    ASSERT_TRUE(rt.outputs[0].npages == 1u);

    page_idx = rt.outputs[0].page_indices[0];
    ASSERT_TRUE(page_idx < rt.npages);
    ASSERT_TRUE(strcmp(rt.pages[page_idx].quantity, "ENERGY") == 0);
    ASSERT_TRUE(rt.pages[page_idx].geometry_idx == 0u);
    ASSERT_TRUE(rt.pages[page_idx].output_idx == 0u);
    ASSERT_TRUE(rt.pages[page_idx].len == 10u);
    ASSERT_TRUE(rt.pages[page_idx].data != NULL);
    ASSERT_TRUE(rt.pages[page_idx].nfilters == 0u);

    ASSERT_TRUE(strcmp(rt.outputs[1].filename, "NB_msh_fluence.bdo") == 0);
    ASSERT_TRUE(strcmp(rt.outputs[1].fileformat, "BDO") == 0);
    ASSERT_TRUE(rt.outputs[1].geometry_idx == 0u);
    ASSERT_TRUE(rt.outputs[1].npages == 3u);

    page_idx = rt.outputs[1].page_indices[2];
    ASSERT_TRUE(page_idx < rt.npages);
    ASSERT_TRUE(strcmp(rt.pages[page_idx].quantity, "FLUENCE") == 0);
    ASSERT_TRUE(rt.pages[page_idx].output_idx == 1u);
    ASSERT_TRUE(rt.pages[page_idx].geometry_idx == 0u);
    ASSERT_TRUE(rt.pages[page_idx].len == 10u);
    ASSERT_TRUE(rt.pages[page_idx].nfilters == 2u);
    ASSERT_TRUE(rt.pages[page_idx].filters[0].filter_idx == 0u);
    ASSERT_TRUE(rt.pages[page_idx].filters[1].filter_idx == 1u);
    ASSERT_TRUE(rt.pages[page_idx].data != NULL);
    ASSERT_TRUE(rt.pages[page_idx].data_var == NULL);
    ASSERT_TRUE(rt.pages[page_idx].data2 == NULL);
    ASSERT_TRUE(rt.pages[page_idx].data2_var == NULL);

    osh_transport_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
}

int main(void) {
    test_compile_fixture_test01_detect();
    return 0;
}
