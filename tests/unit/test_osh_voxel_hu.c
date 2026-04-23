#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "gemca/voxel/osh_gemca2_voxel_hu.h"
#include "openshieldhit/material.h"
#include "openshieldhit/status.h"
#include "test_assert.h"

/* _ct_hu[25] = {-1000,-950,-120,-82,-52,-22,8,19,80,120,...,1600}
 * bin index = lower bracket index in that array.
 * HU=-1000 → bin 0 (air)
 * HU=0     → falls in [-22,8) → bin 5
 * HU=1600  → bin 23 (upper limit, last bin)
 */

/* Verify hu2rho is continuous and monotonically increasing across breakpoints.
 * Checks that the Schneider 2000 piecewise fits join without jumps at
 * HU = -98, 14, 23, 100. */
static void test_hu2rho_continuity(void) {
    float rho_lo;
    float rho_hi;

    /* breakpoint at -98: air/lung → adipose-water */
    rho_lo = osh_gemca_voxel_hu2rho(-99, 0);
    rho_hi = osh_gemca_voxel_hu2rho(-97, 0);
    ASSERT_TRUE(rho_hi > rho_lo); /* density increases with HU */

    /* breakpoint at 14: adipose-water → water transition */
    rho_lo = osh_gemca_voxel_hu2rho(13, 0);
    rho_hi = osh_gemca_voxel_hu2rho(15, 0);
    ASSERT_TRUE(rho_hi > rho_lo);

    /* air density at HU = -1000 */
    ASSERT_TRUE(osh_gemca_voxel_hu2rho(-1000, 0) < 0.002f); /* ~0.00121 g/cm³ */

    /* water at HU = 0 */
    float rho_water = osh_gemca_voxel_hu2rho(0, 0);
    ASSERT_TRUE(rho_water > 0.99f && rho_water < 1.05f);

    /* bone at HU = 1000 */
    ASSERT_TRUE(osh_gemca_voxel_hu2rho(1000, 0) > 1.5f);
}

static void test_hu_lut_boundary_values(void) {
    uint8_t lut[2601];

    osh_gemca_voxel_build_hu_lut(lut);

    ASSERT_TRUE(lut[(-1000) + 1000] == 0); /* air, lowest HU */
    ASSERT_TRUE(lut[(-949) + 1000] == 1);  /* first HU strictly inside second bin */
    ASSERT_TRUE(lut[0 + 1000] == 5);       /* 0 HU: falls in [-22, 8) */
    ASSERT_TRUE(lut[1600 + 1000] == 23);   /* dense bone, upper limit */
}

static void test_hu_lut_matches_hu2idx(void) {
    uint8_t lut[2601];
    int hu;

    osh_gemca_voxel_build_hu_lut(lut);

    hu = -1000;
    while (hu <= 1600) {
        ASSERT_TRUE((int) lut[hu + 1000] == osh_gemca_voxel_hu2idx((int16_t) hu));
        hu++;
    }
}

static void test_schneider_registration_count(void) {
    struct osh_material_workspace *wm = NULL;
    enum osh_status rc;

    rc = osh_material_workspace_create(&wm);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wm != NULL);

    rc = osh_gemca_voxel_register_schneider_materials(wm);
    ASSERT_TRUE(rc == OSH_OK);

    /* 2 reserved (blackhole, vacuum) + 24 Schneider bins */
    ASSERT_TRUE(wm->nmaterials == 26u);

    osh_material_workspace_free(wm);
}

static void test_schneider_registration_names(void) {
    struct osh_material_workspace *wm = NULL;
    enum osh_status rc;

    rc = osh_material_workspace_create(&wm);
    ASSERT_TRUE(rc == OSH_OK);

    rc = osh_gemca_voxel_register_schneider_materials(wm);
    ASSERT_TRUE(rc == OSH_OK);

    ASSERT_TRUE(osh_material_by_name(wm, "schneider_00") != NULL);
    ASSERT_TRUE(osh_material_by_name(wm, "schneider_23") != NULL);
    ASSERT_TRUE(osh_material_by_name(wm, "schneider_24") == NULL);

    osh_material_workspace_free(wm);
}

static void test_schneider_registration_properties(void) {
    struct osh_material_workspace *wm = NULL;
    struct osh_material const *air;
    struct osh_material const *bone;
    enum osh_status rc;

    rc = osh_material_workspace_create(&wm);
    ASSERT_TRUE(rc == OSH_OK);

    rc = osh_gemca_voxel_register_schneider_materials(wm);
    ASSERT_TRUE(rc == OSH_OK);

    air = osh_material_by_name(wm, "schneider_00");
    ASSERT_TRUE(air != NULL);
    ASSERT_TRUE(air->state == OSH_MATERIAL_STATE_GAS);
    ASSERT_TRUE(air->rho > 0.001 && air->rho < 0.002); /* ~0.00121 g/cm^3 */
    ASSERT_TRUE(air->nelements > 0u);

    bone = osh_material_by_name(wm, "schneider_23");
    ASSERT_TRUE(bone != NULL);
    ASSERT_TRUE(bone->state == OSH_MATERIAL_STATE_CONDENSED);
    ASSERT_TRUE(bone->rho > 1.9 && bone->rho < 2.0); /* ~1.93460 g/cm^3 */

    osh_material_workspace_free(wm);
}

static void test_schneider_passes_prepare(void) {
    struct osh_material_workspace *wm = NULL;
    enum osh_status rc;

    rc = osh_material_workspace_create(&wm);
    ASSERT_TRUE(rc == OSH_OK);

    rc = osh_gemca_voxel_register_schneider_materials(wm);
    ASSERT_TRUE(rc == OSH_OK);

    rc = osh_material_workspace_prepare(wm, NULL);
    ASSERT_TRUE(rc == OSH_OK);

    osh_material_workspace_free(wm);
}

static void test_permatassari_hu_lut_boundary_values(void) {
    uint8_t lut[2601];

    osh_gemca_voxel_build_hu_lut_permatassari2020(lut);

    ASSERT_TRUE(lut[(-1000) + 1000] == 0);
    ASSERT_TRUE(lut[0 + 1000] == 14);
    ASSERT_TRUE(lut[1600 + 1000] == 39);
}

static void test_permatassari_hu2rho_reference_point(void) {
    float rho;

    rho = osh_gemca_voxel_hu2rho_permatassari2020(0, 14);
    ASSERT_TRUE(rho > 0.997f && rho < 0.999f); /* expected ~0.998 */

    /* The first bin can go near zero in TOPAS; our implementation clamps it. */
    ASSERT_TRUE(osh_gemca_voxel_hu2rho_permatassari2020(-1000, 0) >= 0.0001f);
}

static void test_rho_lut_schneider_values(void) {
    float lut[2601];

    osh_gemca_voxel_build_rho_lut_schneider2000(lut);

    /* air at HU=-1000: Eq.20 gives 1.03091 + 0.0010297*(-1000) ≈ 0.00121 */
    ASSERT_TRUE(lut[(-1000) + 1000] < 0.002f);
    ASSERT_TRUE(lut[(-1000) + 1000] > 0.0f);

    /* water at HU=0: Eq.21 gives 1.018 + 0.000893*0 = 1.018 */
    ASSERT_TRUE(lut[0 + 1000] > 0.99f && lut[0 + 1000] < 1.05f);

    /* dense bone at HU=1600: Eq.24 gives 1.017 + 0.000592*1600 ≈ 1.964 */
    ASSERT_TRUE(lut[1600 + 1000] > 1.5f);

    /* LUT matches the scalar function at every point */
    {
        int hu;
        for (hu = -1000; hu <= 1600; hu++) {
            ASSERT_TRUE(lut[hu + 1000] == osh_gemca_voxel_hu2rho((int16_t) hu, 0));
        }
    }
}

static void test_rho_lut_permatassari_values(void) {
    float lut[2601];

    osh_gemca_voxel_build_rho_lut_permatassari2020(lut);

    /* water at HU=0, bin 14: factor[14]*(1000+0) ≈ 0.998 */
    ASSERT_TRUE(lut[0 + 1000] > 0.997f && lut[0 + 1000] < 0.999f);

    /* air at HU=-1000 is clamped to >= 0.0001 */
    ASSERT_TRUE(lut[(-1000) + 1000] >= 0.0001f);

    /* dense bone at HU=1600 should have density > 1.5 */
    ASSERT_TRUE(lut[1600 + 1000] > 1.5f);
}

static void test_permatassari_registration_count(void) {
    struct osh_material_workspace *wm = NULL;
    enum osh_status rc;

    rc = osh_material_workspace_create(&wm);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wm != NULL);

    rc = osh_gemca_voxel_register_permatassari_materials(wm);
    ASSERT_TRUE(rc == OSH_OK);

    /* 2 reserved (blackhole, vacuum) + 40 Permatassari bins */
    ASSERT_TRUE(wm->nmaterials == 42u);

    osh_material_workspace_free(wm);
}

static void test_permatassari_registration_names_and_properties(void) {
    struct osh_material_workspace *wm = NULL;
    struct osh_material const *air;
    struct osh_material const *bone;
    enum osh_status rc;

    rc = osh_material_workspace_create(&wm);
    ASSERT_TRUE(rc == OSH_OK);

    rc = osh_gemca_voxel_register_permatassari_materials(wm);
    ASSERT_TRUE(rc == OSH_OK);

    ASSERT_TRUE(osh_material_by_name(wm, "permatassari2020_00") != NULL);
    ASSERT_TRUE(osh_material_by_name(wm, "permatassari2020_39") != NULL);
    ASSERT_TRUE(osh_material_by_name(wm, "permatassari2020_40") == NULL);

    air = osh_material_by_name(wm, "permatassari2020_00");
    ASSERT_TRUE(air != NULL);
    ASSERT_TRUE(air->state == OSH_MATERIAL_STATE_GAS);
    ASSERT_TRUE(air->mean_excitation_energy > 0.0);
    ASSERT_TRUE(air->nelements > 0u);

    bone = osh_material_by_name(wm, "permatassari2020_39");
    ASSERT_TRUE(bone != NULL);
    ASSERT_TRUE(bone->state == OSH_MATERIAL_STATE_CONDENSED);
    ASSERT_TRUE(bone->mean_excitation_energy > 100.0);
    ASSERT_TRUE(bone->nelements > 0u);

    osh_material_workspace_free(wm);
}

int main(void) {
    test_hu2rho_continuity();
    test_hu_lut_boundary_values();
    test_hu_lut_matches_hu2idx();
    test_schneider_registration_count();
    test_schneider_registration_names();
    test_schneider_registration_properties();
    test_schneider_passes_prepare();
    test_permatassari_hu_lut_boundary_values();
    test_permatassari_hu2rho_reference_point();
    test_permatassari_registration_count();
    test_permatassari_registration_names_and_properties();
    test_rho_lut_schneider_values();
    test_rho_lut_permatassari_values();
    return 0;
}
