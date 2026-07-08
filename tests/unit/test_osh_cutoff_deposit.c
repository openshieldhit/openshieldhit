/*
 * Regression test for issue #279: an ion killed at the transport energy cutoff
 * must deposit its residual kinetic energy locally, not delete it.
 *
 * Setup: 150 MeV protons into the 00_minimal water cylinder (z in [0, 20] cm,
 * range ~15.8 cm so every primary stops inside), with nuclear reactions OFF
 * (NUCRE 0) and multiple scattering OFF (MSCAT 0).  With no hadronic channel
 * and no lateral spread, electromagnetic transport must deposit *all* of each
 * primary's kinetic energy inside the scoring mesh: total scored energy per
 * primary == the beam energy, exactly, by energy conservation.
 *
 * TCUT0 is raised to 2 MeV/u so the leak is large and unambiguous: before the
 * fix each stopping proton's exit energy clamps to exactly the cutoff on its
 * penultimate step and is then deleted, discarding 2 MeV per track (2/150 =
 * 1.33% of the beam).  The fix deposits that residual, restoring balance.
 */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "apps/osh/osh_app_osh.h"
#include "openshieldhit/geometry.h"
#include "openshieldhit/material.h"
#include "openshieldhit/scoring.h"
#include "openshieldhit/simulation.h"
#include "openshieldhit/status.h"
#include "test_assert.h"

static int tmp_counter = 0;

static void write_temp_file(char *path_out, size_t cap, char const *content) {
    FILE *fp;
    snprintf(path_out, cap, "osh_cutoff_deposit_%d.tmp", tmp_counter++);
    fp = fopen(path_out, "w");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fputs(content, fp) >= 0);
    fclose(fp);
}

/* Sum the last whitespace-separated column of every data (non-'#') line. */
static double sum_last_column(char const *path) {
    FILE *fp;
    char line[4096];
    double sum;

    sum = 0.0;
    fp = fopen(path, "r");
    ASSERT_TRUE(fp != NULL);
    while (fgets(line, sizeof(line), fp)) {
        char *p;
        char *tok;
        char *last;

        p = line;
        last = NULL;
        while (*p == ' ' || *p == '\t') {
            ++p;
        }
        if (*p == '#' || *p == '\n' || *p == '\r' || *p == '\0') {
            continue;
        }
        for (tok = strtok(p, " \t\r\n"); tok != NULL; tok = strtok(NULL, " \t\r\n")) {
            last = tok;
        }
        if (last) {
            sum += strtod(last, NULL);
        }
    }
    fclose(fp);
    return sum;
}

static void test_cutoff_residual_is_deposited(void) {
    char geo_path[512];
    char mat_path[512];
    char beam_path[512];
    char scoring_path[512];
    char scoring_text[1024];
    char const beam_text[] = "RNDSEED   12345\n"
                             "PRIMARY   1 1\n"
                             "TMAX0     150.0 0.0\n"
                             "BEAMSIGMA 0.0 0.0\n"
                             "BEAMPOS   0.0 0.0 -5.0\n"
                             "NSTAT     2000 -1\n"
                             "DELTAE    0.005\n"
                             "DEMIN     0.025\n"
                             "TCUT0     2.0 150.0\n" /* raised cutoff: 2 MeV/u leak per track */
                             "STRAGG    1\n"         /* exercises the post-straggling exit clamp */
                             "MSCAT     0\n"         /* no lateral spread: all energy stays in the mesh */
                             "NUCRE     0\n";        /* EM only: no hadronic escape channel */
    char const *out_path = "osh_cutoff_deposit_out.dat";
    struct osh_geometry_workspace *geo = NULL;
    struct osh_beam_workspace *beam = NULL;
    struct osh_material_workspace *mat = NULL;
    struct osh_scoring_workspace *scoring = NULL;
    struct osh_simulation *sim = NULL;
    struct osh_results const *results = NULL;
    double energy_per_primary;
    unsigned long long completed;

    snprintf(geo_path, sizeof(geo_path), "%s/tests/cases/00_minimal/geo.dat", OSH_PROJECT_SOURCE_DIR);
    snprintf(mat_path, sizeof(mat_path), "%s/tests/cases/00_minimal/mat.dat", OSH_PROJECT_SOURCE_DIR);
    write_temp_file(beam_path, sizeof(beam_path), beam_text);

    /* One big bin spanning the whole cylinder: the single data line's value is
     * the total deposited energy divided by nstat (a per-primary mean). */
    snprintf(scoring_text,
             sizeof(scoring_text),
             "Geometry Mesh\n"
             "  Name Ebal\n"
             "  X -10.0 10.0 1\n"
             "  Y -10.0 10.0 1\n"
             "  Z 0.0 20.0 1\n"
             "Output\n"
             "  Filename %s\n"
             "  Fileformat TEXT\n"
             "  Geo Ebal\n"
             "  Quantity Energy\n",
             out_path);
    write_temp_file(scoring_path, sizeof(scoring_path), scoring_text);

    ASSERT_TRUE(osh_geometry_setup_from_path(geo_path, NULL, &geo) == OSH_OK);
    ASSERT_TRUE(osh_beam_setup_from_path(beam_path, NULL, &beam) == OSH_OK);
    ASSERT_TRUE(osh_material_setup_from_path(mat_path, NULL, &mat) == OSH_OK);
    ASSERT_TRUE(osh_scoring_setup_from_path(scoring_path, NULL, &scoring) == OSH_OK);

    ASSERT_TRUE(osh_simulation_create(beam, geo, mat, scoring, NULL, &sim) == OSH_OK);
    ASSERT_TRUE(osh_simulation_run(sim) == OSH_OK);
    ASSERT_TRUE(osh_simulation_get_results(sim, &results) == OSH_OK);
    completed = osh_results_completed_nstat(results);
    ASSERT_TRUE(completed > 0ull);

    ASSERT_TRUE(osh_simulation_save(sim) == OSH_OK);
    energy_per_primary = sum_last_column(out_path);

    /* Energy conservation: all 150 MeV of every stopping proton is deposited.
     * The bug leaves this at ~148 MeV (2 MeV cutoff deleted per track); the
     * window is far tighter than that 2 MeV gap. */
    ASSERT_TRUE(fabs(energy_per_primary - 150.0) < 0.05);

    ASSERT_TRUE(osh_simulation_free(sim) == OSH_OK);
    osh_geometry_workspace_free(geo);
    osh_beam_workspace_free(beam);
    osh_material_workspace_free(mat);
    osh_scoring_workspace_free(scoring);
    remove(out_path);
    remove(beam_path);
    remove(scoring_path);
}

int main(void) {
    test_cutoff_residual_is_deposited();
    printf("test_osh_cutoff_deposit: all tests passed\n");
    return 0;
}
