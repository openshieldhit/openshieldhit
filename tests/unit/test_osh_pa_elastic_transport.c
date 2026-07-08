/*
 * Regression test for issue #275: a p+A elastic scatter must update the
 * surviving primary's direction and energy in the particle pool, not only in
 * the scored step.  With the bug, primaries continue straight at undiminished
 * energy after every p+A elastic event, so a pencil beam in an elastic-only
 * run (NUCRE 2, MSCAT off) never leaves a narrow scoring column and its
 * primary fluence at depth stays at the pp-elastic-only level (~0.98 at
 * two-thirds range).  With p+A elastic actually deflecting (~9 deg median lab
 * angle at these energies), the column fluence drops well below that.
 *
 * The bounds are statistical, not bitwise (RNG-stream shifts must not break
 * this test): with 5000 primaries the binomial sigma is ~0.003.  Landmarks
 * for phi_depth (the z = 10..11 cm primary column fluence): ~0.985 with the
 * transport bug (pp-only removal), ~0.96 with the calibrated #277 cross
 * sections (LA150 sigma_R + energy-dependent sigma_el/sigma_R ratio), and
 * ~0.87 with the pre-#277 black-disk sigma_el = sigma_R.  The asserted
 * window (0.90, 0.972) is >4 sigma from both failure modes.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
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
    snprintf(path_out, cap, "osh_pa_elastic_%d.tmp", tmp_counter++);
    fp = fopen(path_out, "w");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fputs(content, fp) >= 0);
    fclose(fp);
}

/* Read the scored value of the mesh bin whose z-column matches @p z_center
 * from a TEXT output ("x y z value" data lines, '#' comments skipped). */
static double read_bin_value(char const *path, double z_center) {
    FILE *fp;
    char line[512];
    double x;
    double y;
    double z;
    double v;
    double found;
    int have;

    fp = fopen(path, "r");
    ASSERT_TRUE(fp != NULL);
    found = 0.0;
    have = 0;
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#') {
            continue;
        }
        if (sscanf(line, "%lf %lf %lf %lf", &x, &y, &z, &v) != 4) {
            continue;
        }
        if (fabs(z - z_center) < 0.25) {
            found = v;
            have = 1;
        }
    }
    fclose(fp);
    ASSERT_TRUE(have == 1);
    return found;
}

/* Elastic-only (NUCRE 2) 150 MeV pencil beam in water, MSCAT and straggling
 * off, primary-proton fluence scored on a 1 cm^2 column: p+A elastic is the
 * only mechanism (beyond the small pp channel) that can move primaries off
 * the axis, so the column fluence at depth measures whether the channel
 * deflects primaries at all. */
static void test_pa_elastic_depletes_narrow_column(void) {
    char geo_path[512];
    char mat_path[512];
    char beam_path[512];
    char scoring_path[512];
    char scoring_text[1024];
    char const *out_path = "osh_pa_elastic_out.dat";
    struct osh_geometry_workspace *geo = NULL;
    struct osh_beam_workspace *beam = NULL;
    struct osh_material_workspace *mat = NULL;
    struct osh_scoring_workspace *scoring = NULL;
    struct osh_simulation *sim = NULL;
    double phi_entrance;
    double phi_depth;

    /* Water cylinder R=10 cm, L=20 cm: laterally wide, so deflected primaries
     * stay in the phantom and only leave the scored column. */
    snprintf(geo_path, sizeof(geo_path), "%s/tests/cases/13_minimal_nucre/geo.dat", OSH_PROJECT_SOURCE_DIR);
    snprintf(mat_path, sizeof(mat_path), "%s/tests/cases/13_minimal_nucre/mat.dat", OSH_PROJECT_SOURCE_DIR);

    write_temp_file(beam_path,
                    sizeof(beam_path),
                    "RNDSEED     20260708           # fixed seed\n"
                    "PRIMARY     1     1            # protons\n"
                    "TMAX0       150.0 0.0          # 150 MeV, range ~15.8 cm in water\n"
                    "BEAMSIGMA   0.0   0.0          # pencil beam\n"
                    "BEAMPOS     0.0   0.0  -5.0\n"
                    "NSTAT       5000  -1\n"
                    "DELTAE      0.005\n"
                    "DEMIN       0.025\n"
                    "STRAGG      0                  # no straggling\n"
                    "MSCAT       0                  # no MCS: nuclear elastic is the only deflection\n"
                    "NUCRE       2                  # elastic only (pp + p+A), no inelastic\n");

    snprintf(scoring_text,
             sizeof(scoring_text),
             "Filter\n"
             "    Name Primary\n"
             "    GEN = 0\n"
             "    Z = 1\n"
             "    A = 1\n"
             "Geometry Mesh\n"
             "    Name Column\n"
             "    X -0.5 0.5 1\n"
             "    Y -0.5 0.5 1\n"
             "    Z  0.0 16.0 16\n"
             "Output\n"
             "    Filename %s\n"
             "    Fileformat TEXT\n"
             "    Geo Column\n"
             "    Quantity Fluence Primary\n",
             out_path);
    write_temp_file(scoring_path, sizeof(scoring_path), scoring_text);

    ASSERT_TRUE(osh_geometry_setup_from_path(geo_path, NULL, &geo) == OSH_OK);
    ASSERT_TRUE(osh_beam_setup_from_path(beam_path, NULL, &beam) == OSH_OK);
    ASSERT_TRUE(osh_material_setup_from_path(mat_path, NULL, &mat) == OSH_OK);
    ASSERT_TRUE(osh_scoring_setup_from_path(scoring_path, NULL, &scoring) == OSH_OK);

    ASSERT_TRUE(osh_simulation_create(beam, geo, mat, scoring, NULL, &sim) == OSH_OK);
    ASSERT_TRUE(osh_simulation_run(sim) == OSH_OK);
    ASSERT_TRUE(osh_simulation_save(sim) == OSH_OK);

    /* Entrance bin (z = 0..1 cm): essentially the full beam.  Two-thirds range
     * (z = 10..11 cm): pp elastic alone leaves ~0.985 of the beam in the column
     * (the #275 bug); working p+A elastic with the calibrated ratio lands near
     * 0.96.  Parse each bin once and print both so a statistical failure shows
     * the measured values. */
    phi_entrance = read_bin_value(out_path, 0.5);
    phi_depth = read_bin_value(out_path, 10.5);
    printf("pa_elastic: phi_entrance=%.4f phi_depth=%.4f\n", phi_entrance, phi_depth);

    ASSERT_TRUE(phi_entrance > 0.97);
    ASSERT_TRUE(phi_entrance < 1.03);

    /* The floor also trips if the elastic magnitude regresses to the black-disk
     * scale (~0.87) or beyond. */
    ASSERT_TRUE(phi_depth < 0.972);
    ASSERT_TRUE(phi_depth > 0.90);

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
    test_pa_elastic_depletes_narrow_column();
    printf("test_osh_pa_elastic_transport: all tests passed\n");
    return 0;
}
