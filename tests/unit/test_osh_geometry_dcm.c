#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_geometry_parse.h"
#include "common/osh_file.h"
#include "common/osh_voxel_order.h"
#include "gemca/osh_gemca2.h"
#include "gemca/runtime/osh_gemca_runtime.h"
#include "openshieldhit/dicom.h"
#include "openshieldhit/geometry.h"
#include "openshieldhit/geometry_defs.h"
#include "openshieldhit/status.h"
#include "openshieldhit/voxel.h"
#include "test_assert.h"

#define CT_DIR OSH_TEST_FIXTURES_DIR "/dicom/DCPT_headphantom"

static void test_dcm_card_populates_vox_body_arguments(void);
static void test_dcm_prepare_compile_propagates_ct_grid(void);
static void test_dcm_transform_matrix_is_orthonormal(void);
static void test_legacy_vox_card_reports_todo(void);
static void assert_rotation_orthonormal(double const t[16], double eps);
static void write_basic_dcm_geo(char const *geo_path);
static void write_dcm_geo_with_angles(char const *geo_path, double gantry_deg, double couch_deg);

static int nearly_equal(double a, double b, double eps) {
    return fabs(a - b) <= eps;
}

int main(void) {
    test_dcm_card_populates_vox_body_arguments();
    test_dcm_prepare_compile_propagates_ct_grid();
    test_dcm_transform_matrix_is_orthonormal();
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
    ASSERT_TRUE(ws->bodies[0].hu != NULL);
    {
        size_t expected_n_hu = (size_t) ct.cols * (size_t) ct.rows * (size_t) ct.n_slices;
        if (OSH_VOXEL_LAYOUT_DEFAULT == OSH_VOXEL_ORDER_MORTON8) {
            size_t Tx = ((size_t) ct.cols + 7u) >> 3u;
            size_t Ty = ((size_t) ct.rows + 7u) >> 3u;
            size_t Tz = ((size_t) ct.n_slices + 7u) >> 3u;
            expected_n_hu = Tx * Ty * Tz * 512u;
        }
        ASSERT_TRUE(ws->bodies[0].n_hu == expected_n_hu);
    }
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
    ASSERT_TRUE(cold_body->hu != NULL);
    ASSERT_TRUE(cold_body->hu == ws->bodies[0].hu);
    ASSERT_TRUE(nearly_equal(cold_body->ct_grid.origin[0], ws->bodies[0].a[0], 1.0e-9));
    ASSERT_TRUE(nearly_equal(cold_body->ct_grid.origin[1], ws->bodies[0].a[1], 1.0e-9));
    ASSERT_TRUE(nearly_equal(cold_body->ct_grid.origin[2], ws->bodies[0].a[2], 1.0e-9));
    ASSERT_TRUE(nearly_equal(cold_body->ct_grid.spacing[0], ws->bodies[0].a[3], 1.0e-9));
    ASSERT_TRUE(nearly_equal(cold_body->ct_grid.spacing[1], ws->bodies[0].a[4], 1.0e-9));
    ASSERT_TRUE(nearly_equal(cold_body->ct_grid.spacing[2], ws->bodies[0].a[5], 1.0e-9));
    ASSERT_TRUE(cold_body->ct_grid.n[0] == (size_t) ct.cols);
    ASSERT_TRUE(cold_body->ct_grid.n[1] == (size_t) ct.rows);
    ASSERT_TRUE(cold_body->ct_grid.n[2] == (size_t) ct.n_slices);

    rc = osh_gemca_compile(ws->prepared, OSH_HU_TABLE_SCHNEIDER, NULL, &rt);
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
    ASSERT_TRUE(rt_body->hu != NULL);
    ASSERT_TRUE(rt_body->hu == cold_body->hu);
    ASSERT_TRUE(rt_body->ct_grid.tile_order == OSH_VOXEL_LAYOUT_DEFAULT);
    ASSERT_TRUE(rt_body->hu[0] == ct.pixels[0]);
    {
        size_t ri = 1u + (size_t) ct.cols * (2u + (size_t) ct.rows * 3u);
        size_t hi = ri;
        if (OSH_VOXEL_LAYOUT_DEFAULT == OSH_VOXEL_ORDER_MORTON8) {
            size_t Tx = ((size_t) ct.cols + 7u) >> 3u;
            size_t Ty = ((size_t) ct.rows + 7u) >> 3u;
            hi = osh_voxel_tile_idx(1u, 2u, 3u, Tx, Ty);
        }
        ASSERT_TRUE(rt_body->hu[hi] == ct.pixels[ri]);
    }

    osh_gemca_runtime_free(&rt);
    osh_dicom_ct_free(&ct);
    osh_fclose(geo);
    (void) osh_geometry_workspace_free(ws);
    (void) remove(geo_path);
}

static void test_dcm_transform_matrix_is_orthonormal(void) {
    double const angle_pairs[][2] = {
        {0.0, 0.0},
        {30.0, -15.0},
        {90.0, 0.0},
        {-45.0, 20.0},
        {180.0, -90.0},
    };
    size_t i;

    for (i = 0u; i < sizeof(angle_pairs) / sizeof(angle_pairs[0]); ++i) {
        char geo_path[64];
        struct oshfile *geo = NULL;
        struct osh_geometry_workspace *ws = NULL;
        struct body const *cold_body;
        enum osh_status rc;

        (void) snprintf(geo_path, sizeof(geo_path), "osh_test_geo_dcm_rot_%u.dat", (unsigned int) i);
        write_dcm_geo_with_angles(geo_path, angle_pairs[i][0], angle_pairs[i][1]);

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
        assert_rotation_orthonormal(cold_body->t, 1.0e-9);

        osh_fclose(geo);
        (void) osh_geometry_workspace_free(ws);
        (void) remove(geo_path);
    }
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
    write_dcm_geo_with_angles(geo_path, 30.0, -15.0);
}

static void write_dcm_geo_with_angles(char const *geo_path, double gantry_deg, double couch_deg) {
    FILE *fp = fopen(geo_path, "w");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fprintf(fp, " 0 0 test\n") > 0);
    ASSERT_TRUE(fprintf(fp, " DCM CTBOX %s %.12g %.12g 1.25 -2.5 3.75\n", CT_DIR, gantry_deg, couch_deg) > 0);
    ASSERT_TRUE(fprintf(fp, " END\n") > 0);
    ASSERT_TRUE(fprintf(fp, " Z001 +CTBOX\n") > 0);
    ASSERT_TRUE(fprintf(fp, " END\n") > 0);
    ASSERT_TRUE(fprintf(fp, " ASSIGNMAT Water Z001\n") > 0);
    fclose(fp);
}

static void assert_rotation_orthonormal(double const t[16], double eps) {
    double const r00 = t[0];
    double const r01 = t[1];
    double const r02 = t[2];
    double const r10 = t[4];
    double const r11 = t[5];
    double const r12 = t[6];
    double const r20 = t[8];
    double const r21 = t[9];
    double const r22 = t[10];
    double const row0_norm = sqrt(r00 * r00 + r01 * r01 + r02 * r02);
    double const row1_norm = sqrt(r10 * r10 + r11 * r11 + r12 * r12);
    double const row2_norm = sqrt(r20 * r20 + r21 * r21 + r22 * r22);
    double const row01_dot = r00 * r10 + r01 * r11 + r02 * r12;
    double const row02_dot = r00 * r20 + r01 * r21 + r02 * r22;
    double const row12_dot = r10 * r20 + r11 * r21 + r12 * r22;
    double const det = r00 * (r11 * r22 - r12 * r21) - r01 * (r10 * r22 - r12 * r20) + r02 * (r10 * r21 - r11 * r20);

    ASSERT_TRUE(nearly_equal(row0_norm, 1.0, eps));
    ASSERT_TRUE(nearly_equal(row1_norm, 1.0, eps));
    ASSERT_TRUE(nearly_equal(row2_norm, 1.0, eps));
    ASSERT_TRUE(nearly_equal(row01_dot, 0.0, eps));
    ASSERT_TRUE(nearly_equal(row02_dot, 0.0, eps));
    ASSERT_TRUE(nearly_equal(row12_dot, 0.0, eps));
    ASSERT_TRUE(nearly_equal(det, 1.0, eps));
}
