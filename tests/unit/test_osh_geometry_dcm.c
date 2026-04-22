#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_geometry_parse.h"
#include "common/osh_file.h"
#include "gemca/osh_gemca2.h"
#include "gemca/runtime/osh_gemca_runtime.h"
#include "openshieldhit/dicom.h"
#include "openshieldhit/geometry.h"
#include "openshieldhit/geometry_defs.h"
#include "openshieldhit/status.h"
#include "test_assert.h"

#define CT_DIR OSH_TEST_FIXTURES_DIR "/dicom/DCPT_headphantom"

static void test_dcm_card_populates_vox_body_arguments(void);
static void test_dcm_prepare_compile_propagates_ct_grid(void);
static void test_legacy_vox_card_reports_todo(void);
static void write_basic_dcm_geo(char const *geo_path);

static int nearly_equal(double a, double b, double eps) {
    return fabs(a - b) <= eps;
}

int main(void) {
    test_dcm_card_populates_vox_body_arguments();
    test_dcm_prepare_compile_propagates_ct_grid();
    test_legacy_vox_card_reports_todo();
    return 0;
}

static void test_dcm_card_populates_vox_body_arguments(void) {
    char const *geo_path = "osh_test_geo_dcm.dat";
    struct oshfile *geo = NULL;
    struct osh_geometry_workspace *ws = NULL;
    struct osh_dicom_ct ct;
    enum osh_status rc;
    write_basic_dcm_geo(geo_path);

    rc = osh_dicom_ct_read(CT_DIR, &ct, NULL);
    ASSERT_TRUE(rc == OSH_OK);

    rc = osh_geometry_workspace_create(&ws);
    ASSERT_TRUE(rc == OSH_OK);
    geo = osh_fopen(geo_path);
    ASSERT_TRUE(geo != NULL);
    rc = osh_geometry_parse(geo, NULL, ws);
    ASSERT_TRUE(rc == OSH_OK);

    ASSERT_TRUE(ws->nbodies == 1u);
    ASSERT_TRUE(ws->bodies[0].type == OSH_GEOMETRY_BODY_VOX);
    ASSERT_TRUE(ws->bodies[0].name != NULL);
    ASSERT_TRUE(strcmp(ws->bodies[0].name, "CTBOX") == 0);
    ASSERT_TRUE(ws->bodies[0].a != NULL);
    ASSERT_TRUE(ws->bodies[0].na == 14);

    ASSERT_TRUE(nearly_equal(ws->bodies[0].a[0], 0.0, 1.0e-9));
    ASSERT_TRUE(nearly_equal(ws->bodies[0].a[1], 0.0, 1.0e-9));
    ASSERT_TRUE(nearly_equal(ws->bodies[0].a[2], 0.0, 1.0e-9));
    ASSERT_TRUE(nearly_equal(ws->bodies[0].a[3], 0.1 * ct.pixel_spacing[1], 1.0e-6));
    ASSERT_TRUE(nearly_equal(ws->bodies[0].a[4], 0.1 * ct.pixel_spacing[0], 1.0e-6));
    ASSERT_TRUE(nearly_equal(ws->bodies[0].a[5], 0.1 * ct.slice_spacing, 1.0e-6));
    ASSERT_TRUE(nearly_equal(ws->bodies[0].a[6], (double) ct.cols, 1.0e-9));
    ASSERT_TRUE(nearly_equal(ws->bodies[0].a[7], (double) ct.rows, 1.0e-9));
    ASSERT_TRUE(nearly_equal(ws->bodies[0].a[8], (double) ct.n_slices, 1.0e-9));
    ASSERT_TRUE(nearly_equal(ws->bodies[0].a[9], 30.0, 1.0e-9));
    ASSERT_TRUE(nearly_equal(ws->bodies[0].a[10], -15.0, 1.0e-9));
    ASSERT_TRUE(nearly_equal(ws->bodies[0].a[11], 1.25 - 0.5 * ws->bodies[0].a[3], 1.0e-9));
    ASSERT_TRUE(nearly_equal(ws->bodies[0].a[12], -2.5 - 0.5 * ws->bodies[0].a[4], 1.0e-9));
    ASSERT_TRUE(nearly_equal(ws->bodies[0].a[13], 3.75 - 0.5 * ws->bodies[0].a[5], 1.0e-9));

    osh_dicom_ct_free(&ct);
    osh_fclose(geo);
    (void) osh_geometry_workspace_free(ws);
    (void) remove(geo_path);
}

static void test_dcm_prepare_compile_propagates_ct_grid(void) {
    char const *geo_path = "osh_test_geo_dcm_compile.dat";
    struct oshfile *geo = NULL;
    struct osh_geometry_workspace *ws = NULL;
    struct osh_gemca_runtime rt;
    struct body const *cold_body;
    struct gemca_rt_body const *rt_body;
    struct osh_dicom_ct ct;
    enum osh_status rc;

    memset(&rt, 0, sizeof(rt));
    write_basic_dcm_geo(geo_path);

    rc = osh_dicom_ct_read(CT_DIR, &ct, NULL);
    ASSERT_TRUE(rc == OSH_OK);

    rc = osh_geometry_workspace_create(&ws);
    ASSERT_TRUE(rc == OSH_OK);
    geo = osh_fopen(geo_path);
    ASSERT_TRUE(geo != NULL);
    rc = osh_geometry_parse(geo, NULL, ws);
    ASSERT_TRUE(rc == OSH_OK);

    rc = osh_geometry_workspace_prepare(ws, NULL);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(ws->prepared != NULL);
    ASSERT_TRUE(ws->prepared->nbodies == 1u);

    cold_body = ws->prepared->bodies[0];
    ASSERT_TRUE(cold_body != NULL);
    ASSERT_TRUE(cold_body->type == OSH_GEOMETRY_BODY_VOX);
    ASSERT_TRUE(nearly_equal(cold_body->ct_grid.origin[0], ws->bodies[0].a[0], 1.0e-9));
    ASSERT_TRUE(nearly_equal(cold_body->ct_grid.origin[1], ws->bodies[0].a[1], 1.0e-9));
    ASSERT_TRUE(nearly_equal(cold_body->ct_grid.origin[2], ws->bodies[0].a[2], 1.0e-9));
    ASSERT_TRUE(nearly_equal(cold_body->ct_grid.spacing[0], ws->bodies[0].a[3], 1.0e-9));
    ASSERT_TRUE(nearly_equal(cold_body->ct_grid.spacing[1], ws->bodies[0].a[4], 1.0e-9));
    ASSERT_TRUE(nearly_equal(cold_body->ct_grid.spacing[2], ws->bodies[0].a[5], 1.0e-9));
    ASSERT_TRUE(cold_body->ct_grid.n[0] == (size_t) ct.cols);
    ASSERT_TRUE(cold_body->ct_grid.n[1] == (size_t) ct.rows);
    ASSERT_TRUE(cold_body->ct_grid.n[2] == (size_t) ct.n_slices);

    rc = osh_gemca_compile(ws->prepared, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(rt.nbodies == 1u);

    rt_body = &rt.bodies[0];
    ASSERT_TRUE(nearly_equal(rt_body->ct_grid.origin[0], cold_body->ct_grid.origin[0], 1.0e-9));
    ASSERT_TRUE(nearly_equal(rt_body->ct_grid.origin[1], cold_body->ct_grid.origin[1], 1.0e-9));
    ASSERT_TRUE(nearly_equal(rt_body->ct_grid.origin[2], cold_body->ct_grid.origin[2], 1.0e-9));
    ASSERT_TRUE(nearly_equal(rt_body->ct_grid.spacing[0], cold_body->ct_grid.spacing[0], 1.0e-9));
    ASSERT_TRUE(nearly_equal(rt_body->ct_grid.spacing[1], cold_body->ct_grid.spacing[1], 1.0e-9));
    ASSERT_TRUE(nearly_equal(rt_body->ct_grid.spacing[2], cold_body->ct_grid.spacing[2], 1.0e-9));
    ASSERT_TRUE(rt_body->ct_grid.n[0] == cold_body->ct_grid.n[0]);
    ASSERT_TRUE(rt_body->ct_grid.n[1] == cold_body->ct_grid.n[1]);
    ASSERT_TRUE(rt_body->ct_grid.n[2] == cold_body->ct_grid.n[2]);
    ASSERT_TRUE(rt_body->hu == NULL);

    osh_gemca_runtime_free(&rt);
    osh_dicom_ct_free(&ct);
    osh_fclose(geo);
    (void) osh_geometry_workspace_free(ws);
    (void) remove(geo_path);
}

static void test_legacy_vox_card_reports_todo(void) {
    char const *geo_path = "osh_test_geo_legacy_vox.dat";
    struct oshfile *geo = NULL;
    struct osh_geometry_workspace *ws = NULL;
    enum osh_status rc;
    FILE *fp = fopen(geo_path, "w");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fprintf(fp, " 0 0 test\n") > 0);
    ASSERT_TRUE(fprintf(fp, " VOX V1 dummy.ctx 0.0 0.0\n") > 0);
    ASSERT_TRUE(fprintf(fp, " END\n") > 0);
    ASSERT_TRUE(fprintf(fp, " Z001 +V1\n") > 0);
    ASSERT_TRUE(fprintf(fp, " END\n") > 0);
    ASSERT_TRUE(fprintf(fp, " ASSIGNMAT Water Z001\n") > 0);
    fclose(fp);

    rc = osh_geometry_workspace_create(&ws);
    ASSERT_TRUE(rc == OSH_OK);
    geo = osh_fopen(geo_path);
    ASSERT_TRUE(geo != NULL);
    rc = osh_geometry_parse(geo, NULL, ws);
    ASSERT_TRUE(rc == OSH_EPARSE);

    osh_fclose(geo);
    (void) osh_geometry_workspace_free(ws);
    (void) remove(geo_path);
}

static void write_basic_dcm_geo(char const *geo_path) {
    FILE *fp = fopen(geo_path, "w");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fprintf(fp, " 0 0 test\n") > 0);
    ASSERT_TRUE(fprintf(fp, " DCM CTBOX %s 30.0 -15.0 1.25 -2.5 3.75\n", CT_DIR) > 0);
    ASSERT_TRUE(fprintf(fp, " END\n") > 0);
    ASSERT_TRUE(fprintf(fp, " Z001 +CTBOX\n") > 0);
    ASSERT_TRUE(fprintf(fp, " END\n") > 0);
    ASSERT_TRUE(fprintf(fp, " ASSIGNMAT Water Z001\n") > 0);
    fclose(fp);
}
