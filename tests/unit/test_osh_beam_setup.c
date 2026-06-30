#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_app_osh.h"
#include "openshieldhit/status.h"
#include "particle/osh_particle_pdg.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

static int _tmp_counter = 0;

static void _write_temp_file(char *path, size_t path_cap, char const *content) {
    FILE *fp;

    snprintf(path, path_cap, "osh_test_%d.tmp", _tmp_counter++);
    fp = fopen(path, "w");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fputs(content, fp) >= 0);
    ASSERT_TRUE(fclose(fp) == 0);
}

static void test_beam_spots_set_replace_and_validate(void) {
    struct osh_beam_workspace *wb = NULL;
    struct osh_beam_spot spots_a[2];
    struct osh_beam_spot spots_b[1];
    struct osh_beam_spot invalid[1];
    int rc;

    memset(spots_a, 0, sizeof(spots_a));
    memset(spots_b, 0, sizeof(spots_b));
    memset(invalid, 0, sizeof(invalid));

    rc = osh_beam_workspace_create(&wb);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wb != NULL);

    spots_a[0].shape = OSH_BEAM_SHAPE_PENCIL;
    spots_a[0].div[0] = 0.001;
    spots_a[0].t0 = 100.0;
    spots_a[0].wt = 1.0;
    spots_a[1].shape = OSH_BEAM_SHAPE_GAUSSIAN;
    spots_a[1].size[0] = 0.2;
    spots_a[1].size[1] = 0.3;
    spots_a[1].t0 = 120.0;
    spots_a[1].wt = 2.0;

    ASSERT_TRUE(osh_beam_spots_set(wb, spots_a, 2u) == OSH_OK);
    ASSERT_TRUE(wb->nspots == 2u);
    ASSERT_TRUE(wb->spots != NULL);
    ASSERT_TRUE(wb->spots[0].spot_id == 1u);
    ASSERT_TRUE(wb->spots[1].spot_id == 2u);
    ASSERT_TRUE(wb->shared.use_div == 1);
    ASSERT_TRUE(fabs(wb->spots[1].size[1] - 0.3) < 1e-12);

    spots_b[0].shape = OSH_BEAM_SHAPE_PENCIL;
    spots_b[0].p[2] = -25.0;
    spots_b[0].t0 = 80.0;
    spots_b[0].wt = 4.0;

    ASSERT_TRUE(osh_beam_spots_set(wb, spots_b, 1u) == OSH_OK);
    ASSERT_TRUE(wb->nspots == 1u);
    ASSERT_TRUE(wb->shared.use_div == 0);
    ASSERT_TRUE(wb->spots[0].spot_id == 1u);
    ASSERT_TRUE(fabs(wb->spots[0].p[2] + 25.0) < 1e-12);

    invalid[0].shape = (char) 99;
    invalid[0].wt = 1.0;
    ASSERT_TRUE(osh_beam_spots_set(wb, invalid, 1u) == OSH_EINVAL);

    ASSERT_TRUE(osh_beam_workspace_free(wb) == OSH_OK);
}

static void test_setup_single_spot_from_beamdat(void) {
    char beam_path[512];
    char beam_text[512];
    struct osh_beam_workspace *wb = NULL;
    int rc;

    snprintf(beam_text,
             sizeof beam_text,
             "PRIMARY proton\n"
             "TMAX0 120.0 1.5\n"
             "BEAMPOS 1.0 2.0 -50.0\n"
             "BEAMSIGMA 0.0 0.0\n"
             "BEAMDIV 2.0 3.0 0.0\n");
    _write_temp_file(beam_path, sizeof(beam_path), beam_text);

    rc = osh_beam_setup_from_path(beam_path, NULL, &wb);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wb != NULL);
    ASSERT_TRUE(wb->nspots == 1);
    ASSERT_TRUE(fabs(wb->spots[0].p[0] - 1.0) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].p[1] - 2.0) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].p[2] + 50.0) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].t0 - 120.0) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].tsigma - 1.5) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].div[0] - 0.002) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].div[1] - 0.003) < 1e-12);

    ASSERT_TRUE(osh_beam_workspace_free(wb) == OSH_OK);
    ASSERT_TRUE(remove(beam_path) == 0);
}

static void test_setup_beamdiv_without_focus_defaults_to_zero(void) {
    char beam_path[512];
    char beam_text[512];
    struct osh_beam_workspace *wb = NULL;
    int rc;

    snprintf(beam_text,
             sizeof beam_text,
             "PRIMARY proton\n"
             "TMAX0 120.0 0.0\n"
             "BEAMPOS 0.0 0.0 -10.0\n"
             "BEAMDIV 2.0 3.0\n");
    _write_temp_file(beam_path, sizeof(beam_path), beam_text);

    rc = osh_beam_setup_from_path(beam_path, NULL, &wb);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wb != NULL);
    ASSERT_TRUE(fabs(wb->spots[0].div[0] - 0.002) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].div[1] - 0.003) < 1e-12);
    ASSERT_TRUE(fabs(wb->shared.focus) < 1e-12);
    ASSERT_TRUE(wb->shared.use_div == 1);

    ASSERT_TRUE(osh_beam_workspace_free(wb) == OSH_OK);
    ASSERT_TRUE(remove(beam_path) == 0);
}

static void test_setup_spotlist_replaces_template_and_inherits_defaults(void) {
    char beam_path[512];
    char spot_path[512];
    char beam_text[1024];
    char spot_text[1024];
    struct osh_beam_workspace *wb = NULL;
    double fwhm_x0;
    double fwhm_y0;
    double fwhm_x1;
    double fwhm_y1;
    int rc;

    fwhm_x0 = 0.2 * (2.0 * sqrt(2.0 * log(2.0)));
    fwhm_y0 = 0.1 * (2.0 * sqrt(2.0 * log(2.0)));
    fwhm_x1 = 0.3 * (2.0 * sqrt(2.0 * log(2.0)));
    fwhm_y1 = 0.4 * (2.0 * sqrt(2.0 * log(2.0)));

    snprintf(spot_text,
             sizeof spot_text,
             "# E_GeV dE_GeV X_cm Y_cm FWHMx_cm FWHMy_cm Weight\n"
             "0.150 0.002 1.25 -2.5 %.17g %.17g 1000\n"
             "0.175 0.003 -4.0 3.5 %.17g %.17g 2000\n",
             fwhm_x0,
             fwhm_y0,
             fwhm_x1,
             fwhm_y1);
    _write_temp_file(spot_path, sizeof(spot_path), spot_text);

    snprintf(beam_text,
             sizeof beam_text,
             "PRIMARY proton\n"
             "TMAX0 200.0 1.5\n"
             "BEAMPOS 9.0 8.0 -35.0\n"
             "BEAMSIGMA 0.0 0.0\n"
             "BEAMDIV 4.0 6.0 0.0\n"
             "USECBEAM %s\n",
             spot_path);
    _write_temp_file(beam_path, sizeof(beam_path), beam_text);

    rc = osh_beam_setup_from_path(beam_path, NULL, &wb);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wb != NULL);
    ASSERT_TRUE(wb->nspots == 2);

    ASSERT_TRUE(fabs(wb->spots[0].p[0] - 1.25) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].p[1] + 2.5) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].p[2] + 35.0) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].t0 - 150.0) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].tsigma - 2.0) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].size[0] - 0.2) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].size[1] - 0.1) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].div[0] - 0.004) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].div[1] - 0.006) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].cor[0]) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].cor[1]) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].wt - 1000.0) < 1e-12);
    ASSERT_TRUE(wb->spots[0].spot_id == 1u);

    ASSERT_TRUE(fabs(wb->spots[1].p[0] + 4.0) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[1].p[1] - 3.5) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[1].p[2] + 35.0) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[1].t0 - 175.0) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[1].tsigma - 3.0) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[1].size[0] - 0.3) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[1].size[1] - 0.4) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[1].div[0] - 0.004) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[1].div[1] - 0.006) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[1].wt - 2000.0) < 1e-12);
    ASSERT_TRUE(wb->spots[1].spot_id == 2u);

    ASSERT_TRUE(wb->shared.use_div == 1);
    ASSERT_TRUE(osh_beam_workspace_free(wb) == OSH_OK);
    ASSERT_TRUE(remove(beam_path) == 0);
    ASSERT_TRUE(remove(spot_path) == 0);
}

static void test_setup_spotlist_overrides_optional_columns(void) {
    char beam_path[512];
    char spot_path[512];
    char beam_text[1024];
    char spot_text[1024];
    struct osh_beam_workspace *wb = NULL;
    double fwhm_x;
    double fwhm_y;
    int rc;

    fwhm_x = 0.25 * (2.0 * sqrt(2.0 * log(2.0)));
    fwhm_y = 0.35 * (2.0 * sqrt(2.0 * log(2.0)));

    snprintf(spot_text, sizeof spot_text, "0.125 0.001 2.0 -1.0 %.17g %.17g 10.0 20.0 0.3 -0.4 42\n", fwhm_x, fwhm_y);
    _write_temp_file(spot_path, sizeof(spot_path), spot_text);

    snprintf(beam_text,
             sizeof beam_text,
             "PRIMARY proton\n"
             "TMAX0 220.0 4.0\n"
             "BEAMPOS 0.0 0.0 -25.0\n"
             "BEAMSIGMA 0.0 0.0\n"
             "BEAMDIV 4.0 6.0 0.0\n"
             "USECBEAM %s\n",
             spot_path);
    _write_temp_file(beam_path, sizeof(beam_path), beam_text);

    rc = osh_beam_setup_from_path(beam_path, NULL, &wb);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wb != NULL);
    ASSERT_TRUE(wb->nspots == 1);
    ASSERT_TRUE(fabs(wb->spots[0].p[0] - 2.0) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].p[1] + 1.0) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].p[2] + 25.0) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].t0 - 125.0) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].tsigma - 1.0) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].size[0] - 0.25) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].size[1] - 0.35) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].div[0] - 0.01) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].div[1] - 0.02) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].cor[0] - 0.3) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].cor[1] + 0.4) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].wt - 42.0) < 1e-12);

    ASSERT_TRUE(osh_beam_workspace_free(wb) == OSH_OK);
    ASSERT_TRUE(remove(beam_path) == 0);
    ASSERT_TRUE(remove(spot_path) == 0);
}

static void test_setup_primary_name_resolves_particle(void) {
    char beam_path[512];
    char beam_text[512];
    struct osh_beam_workspace *wb = NULL;
    int rc;

    snprintf(beam_text, sizeof beam_text, "PRIMARY proton\nTMAX0 120.0 0.0\nBEAMPOS 0.0 0.0 -10.0\n");
    _write_temp_file(beam_path, sizeof(beam_path), beam_text);

    rc = osh_beam_setup_from_path(beam_path, NULL, &wb);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wb != NULL);
    ASSERT_TRUE(wb->has_primary == 1);
    ASSERT_TRUE(wb->primary.pdg == OSH_PART_PDG_PROTON);
    ASSERT_TRUE(fabs(wb->spots[0].p0) > 0.0);

    ASSERT_TRUE(osh_beam_workspace_free(wb) == OSH_OK);
    ASSERT_TRUE(remove(beam_path) == 0);
}

static void test_setup_primary_pdg_resolves_particle(void) {
    char beam_path[512];
    char beam_text[512];
    struct osh_beam_workspace *wb = NULL;
    int rc;

    snprintf(beam_text, sizeof beam_text, "PRIMARY 2212\nTMAX0 120.0 0.0\nBEAMPOS 0.0 0.0 -10.0\n");
    _write_temp_file(beam_path, sizeof(beam_path), beam_text);

    rc = osh_beam_setup_from_path(beam_path, NULL, &wb);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wb != NULL);
    ASSERT_TRUE(wb->has_primary == 1);
    ASSERT_TRUE(wb->primary.pdg == OSH_PART_PDG_PROTON);
    ASSERT_TRUE(fabs(wb->spots[0].p0) > 0.0);

    ASSERT_TRUE(osh_beam_workspace_free(wb) == OSH_OK);
    ASSERT_TRUE(remove(beam_path) == 0);
}

static void test_setup_primary_za_resolves_ion(void) {
    char beam_path[512];
    char beam_text[512];
    struct osh_beam_workspace *wb = NULL;
    int rc;

    snprintf(beam_text, sizeof beam_text, "PRIMARY 6 12\nTMAX0 400.0 1.0\nBEAMPOS 0.0 0.0 -10.0\n");
    _write_temp_file(beam_path, sizeof(beam_path), beam_text);

    rc = osh_beam_setup_from_path(beam_path, NULL, &wb);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wb != NULL);
    ASSERT_TRUE(wb->has_primary == 1);
    ASSERT_TRUE(wb->primary.pdg == 1000060120);
    ASSERT_TRUE(wb->primary.z == 6);
    ASSERT_TRUE(wb->primary.a == 12);
    ASSERT_TRUE(fabs(wb->spots[0].t0 - 4800.0) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].tsigma - 12.0) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].p0) > 0.0);

    ASSERT_TRUE(osh_beam_workspace_free(wb) == OSH_OK);
    ASSERT_TRUE(remove(beam_path) == 0);
}

static void test_setup_primary_invalid_returns_einval(void) {
    char beam_path[512];
    char beam_text[512];
    struct osh_beam_workspace *wb = NULL;
    int rc;

    snprintf(beam_text, sizeof beam_text, "PRIMARY nosuchparticle\nTMAX0 120.0 0.0\nBEAMPOS 0.0 0.0 -10.0\n");
    _write_temp_file(beam_path, sizeof(beam_path), beam_text);

    rc = osh_beam_setup_from_path(beam_path, NULL, &wb);

    ASSERT_TRUE(rc == OSH_EINVAL);
    ASSERT_TRUE(wb == NULL);

    ASSERT_TRUE(remove(beam_path) == 0);
}

static void test_setup_missing_primary_returns_einval(void) {
    char beam_path[512];
    char beam_text[512];
    struct osh_beam_workspace *wb = NULL;
    int rc;

    snprintf(beam_text, sizeof beam_text, "TMAX0 120.0 0.0\nBEAMPOS 0.0 0.0 -10.0\n");
    _write_temp_file(beam_path, sizeof(beam_path), beam_text);

    rc = osh_beam_setup_from_path(beam_path, NULL, &wb);

    ASSERT_TRUE(rc == OSH_EINVAL);
    ASSERT_TRUE(wb == NULL);

    ASSERT_TRUE(remove(beam_path) == 0);
}

static void test_setup_unknown_key_returns_eparse(void) {
    char beam_path[512];
    char beam_text[512];
    struct osh_beam_workspace *wb = NULL;
    int rc;

    snprintf(beam_text, sizeof beam_text, "PRIMARY proton\nTMAX0 120.0 0.0\nBANANA 1 2 3\nBEAMPOS 0.0 0.0 -10.0\n");
    _write_temp_file(beam_path, sizeof(beam_path), beam_text);

    rc = osh_beam_setup_from_path(beam_path, NULL, &wb);

    ASSERT_TRUE(rc == OSH_EPARSE);
    ASSERT_TRUE(wb == NULL);

    ASSERT_TRUE(remove(beam_path) == 0);
}

static void test_setup_beamsigma_single_value_sets_symmetric_xy(void) {
    char beam_path[512];
    char beam_text[512];
    struct osh_beam_workspace *wb = NULL;
    int rc;

    snprintf(beam_text, sizeof beam_text, "PRIMARY proton\nTMAX0 120.0 0.0\nBEAMPOS 0.0 0.0 -10.0\nBEAMSIGMA 1.0\n");
    _write_temp_file(beam_path, sizeof(beam_path), beam_text);

    rc = osh_beam_setup_from_path(beam_path, NULL, &wb);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wb != NULL);
    ASSERT_TRUE(wb->spots[0].shape == OSH_BEAM_SHAPE_GAUSSIAN);
    ASSERT_TRUE(fabs(wb->spots[0].size[0] - 1.0) < 1e-12);
    ASSERT_TRUE(fabs(wb->spots[0].size[1] - 1.0) < 1e-12);

    ASSERT_TRUE(osh_beam_workspace_free(wb) == OSH_OK);
    ASSERT_TRUE(remove(beam_path) == 0);
}

static void test_setup_nstat_single_value_defaults_nsave_to_zero(void) {
    char beam_path[512];
    char beam_text[512];
    struct osh_beam_workspace *wb = NULL;
    int rc;

    snprintf(beam_text, sizeof beam_text, "PRIMARY proton\nTMAX0 120.0 0.0\nBEAMPOS 0.0 0.0 -10.0\nNSTAT 1234\n");
    _write_temp_file(beam_path, sizeof(beam_path), beam_text);

    rc = osh_beam_setup_from_path(beam_path, NULL, &wb);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wb != NULL);
    ASSERT_TRUE(wb->nstat == 1234u);
    ASSERT_TRUE(wb->nsave == 0u);

    ASSERT_TRUE(osh_beam_workspace_free(wb) == OSH_OK);
    ASSERT_TRUE(remove(beam_path) == 0);
}

static void test_setup_nstat_negative_save_disables_nsave(void) {
    char beam_path[512];
    char beam_text[512];
    struct osh_beam_workspace *wb = NULL;
    int rc;

    snprintf(beam_text, sizeof beam_text, "PRIMARY proton\nTMAX0 120.0 0.0\nBEAMPOS 0.0 0.0 -10.0\nNSTAT 1000 -1\n");
    _write_temp_file(beam_path, sizeof(beam_path), beam_text);

    rc = osh_beam_setup_from_path(beam_path, NULL, &wb);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wb != NULL);
    ASSERT_TRUE(wb->nstat == 1000u);
    ASSERT_TRUE(wb->nsave == 0u);

    ASSERT_TRUE(osh_beam_workspace_free(wb) == OSH_OK);
    ASSERT_TRUE(remove(beam_path) == 0);
}

static void test_setup_nucre_modes_set_independent_flags(void) {
    char beam_path[512];
    char beam_text[512];
    int mode;

    for (mode = 0; mode <= 3; ++mode) {
        struct osh_beam_workspace *wb = NULL;
        int rc;

        snprintf(beam_text,
                 sizeof beam_text,
                 "PRIMARY proton\n"
                 "TMAX0 120.0 0.0\n"
                 "BEAMPOS 0.0 0.0 -10.0\n"
                 "NUCRE %d\n",
                 mode);
        _write_temp_file(beam_path, sizeof(beam_path), beam_text);

        rc = osh_beam_setup_from_path(beam_path, NULL, &wb);

        ASSERT_TRUE(rc == OSH_OK);
        ASSERT_TRUE(wb != NULL);
        ASSERT_TRUE(wb->nuclear_inelastic == (char) (mode == 1 || mode == 3));
        ASSERT_TRUE(wb->nuclear_elastic == (char) (mode == 1 || mode == 2));

        ASSERT_TRUE(osh_beam_workspace_free(wb) == OSH_OK);
        ASSERT_TRUE(remove(beam_path) == 0);
    }
}

static void test_setup_maxtime_sets_wall_budget(void) {
    char beam_path[512];
    char beam_text[512];
    struct osh_beam_workspace *wb = NULL;
    int rc;

    /* MAXTIME accepts the same duration grammar as --max-time; "30m" → 1800 s. */
    snprintf(beam_text,
             sizeof beam_text,
             "PRIMARY proton\nTMAX0 120.0 0.0\nBEAMPOS 0.0 0.0 -10.0\nNSTAT 1000\nMAXTIME 30m\n");
    _write_temp_file(beam_path, sizeof(beam_path), beam_text);

    rc = osh_beam_setup_from_path(beam_path, NULL, &wb);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wb != NULL);
    ASSERT_TRUE(wb->wall_budget_s == 1800.0);

    ASSERT_TRUE(osh_beam_workspace_free(wb) == OSH_OK);
    ASSERT_TRUE(remove(beam_path) == 0);
}

static void test_setup_maxtime_bare_seconds(void) {
    char beam_path[512];
    char beam_text[512];
    struct osh_beam_workspace *wb = NULL;
    int rc;

    snprintf(beam_text,
             sizeof beam_text,
             "PRIMARY proton\nTMAX0 120.0 0.0\nBEAMPOS 0.0 0.0 -10.0\nNSTAT 1000\nMAXTIME 3600\n");
    _write_temp_file(beam_path, sizeof(beam_path), beam_text);

    rc = osh_beam_setup_from_path(beam_path, NULL, &wb);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wb != NULL);
    ASSERT_TRUE(wb->wall_budget_s == 3600.0);

    ASSERT_TRUE(osh_beam_workspace_free(wb) == OSH_OK);
    ASSERT_TRUE(remove(beam_path) == 0);
}

static void test_setup_no_maxtime_defaults_to_unlimited(void) {
    char beam_path[512];
    char beam_text[512];
    struct osh_beam_workspace *wb = NULL;
    int rc;

    snprintf(beam_text, sizeof beam_text, "PRIMARY proton\nTMAX0 120.0 0.0\nBEAMPOS 0.0 0.0 -10.0\nNSTAT 1000\n");
    _write_temp_file(beam_path, sizeof(beam_path), beam_text);

    rc = osh_beam_setup_from_path(beam_path, NULL, &wb);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(wb != NULL);
    ASSERT_TRUE(wb->wall_budget_s == 0.0); /* 0 = unlimited */

    ASSERT_TRUE(osh_beam_workspace_free(wb) == OSH_OK);
    ASSERT_TRUE(remove(beam_path) == 0);
}

static void test_setup_maxtime_invalid_returns_eparse(void) {
    char beam_path[512];
    char beam_text[512];
    struct osh_beam_workspace *wb = NULL;
    int rc;

    snprintf(beam_text,
             sizeof beam_text,
             "PRIMARY proton\nTMAX0 120.0 0.0\nBEAMPOS 0.0 0.0 -10.0\nNSTAT 1000\nMAXTIME notaduration\n");
    _write_temp_file(beam_path, sizeof(beam_path), beam_text);

    rc = osh_beam_setup_from_path(beam_path, NULL, &wb);

    ASSERT_TRUE(rc == OSH_EPARSE);

    ASSERT_TRUE(remove(beam_path) == 0);
}

int main(void) {
    test_beam_spots_set_replace_and_validate();
    test_setup_single_spot_from_beamdat();
    test_setup_beamdiv_without_focus_defaults_to_zero();
    test_setup_spotlist_replaces_template_and_inherits_defaults();
    test_setup_spotlist_overrides_optional_columns();
    test_setup_primary_name_resolves_particle();
    test_setup_primary_pdg_resolves_particle();
    test_setup_primary_za_resolves_ion();
    test_setup_primary_invalid_returns_einval();
    test_setup_missing_primary_returns_einval();
    test_setup_unknown_key_returns_eparse();
    test_setup_beamsigma_single_value_sets_symmetric_xy();
    test_setup_nstat_single_value_defaults_nsave_to_zero();
    test_setup_nstat_negative_save_disables_nsave();
    test_setup_nucre_modes_set_independent_flags();
    test_setup_maxtime_sets_wall_budget();
    test_setup_maxtime_bare_seconds();
    test_setup_no_maxtime_defaults_to_unlimited();
    test_setup_maxtime_invalid_returns_eparse();
    return 0;
}
