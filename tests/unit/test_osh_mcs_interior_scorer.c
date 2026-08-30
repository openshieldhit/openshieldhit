/*
 * Regression test for issue #325: a detect.dat scoring surface placed strictly
 * inside one homogeneous zone must see the multiple scattering accumulated
 * upstream of it.
 *
 * MCS is sampled once per transport substep — at a random hinge, or (for a
 * boundary-limited substep) at the substep end.  A long, low-density zone used
 * to be crossed in a single boundary-limited substep, so the whole zone's
 * deflection was applied only where the particle left the zone: every scorer
 * inside recorded a perfectly unscattered pencil beam.  The step-length cap on
 * the per-step lateral displacement (OSH_TRANSPORT_MCS_LATERAL_MAX_CM) splits
 * such a crossing into hinged substeps, so the spread builds up along the path.
 *
 * Setup: 70.2 MeV protons launched into 2 m of air, with a radial fluence
 * scorer 1.5 m in.  Straggling and nuclear reactions are off and the beam is a
 * true pencil (BEAMSIGMA/BEAMDIV zero), so the scored width is transport-level
 * MCS and nothing else.  MSCAT 1 (Gaussian) is deliberate: the random hinge is
 * exactly additive over substeps for a Gaussian, so the scored width must not
 * depend on how the path happens to be cut into substeps.
 *
 * The same physics is run three ways and must agree:
 *   - one air zone, DELTAE 0.05   (the reported failing configuration)
 *   - one air zone, DELTAE 0.005  (fine steps: the bug hides when DELTAE alone
 *                                  already splits the crossing)
 *   - five air zones, DELTAE 0.05 (the same 2 m subdivided into real zones,
 *                                  which is how the reporter worked around it)
 * Before the fix these give sigma = 0.088 / 0.568 / 0.447 cm; after, all three
 * land on the Highland value ~0.56 cm.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_app_osh.h"
#include "openshieldhit/const.h"
#include "openshieldhit/geometry.h"
#include "openshieldhit/material.h"
#include "openshieldhit/scoring.h"
#include "openshieldhit/simulation.h"
#include "openshieldhit/status.h"
#include "test_assert.h"

/* Scoring plane: a thin slab 150 cm into the air, binned radially. */
#define PLANE_Z_CM 150.0
#define PLANE_NR 20
#define PLANE_R_MAX_CM 5.0

/* Highland projected width for this deck (70.2 MeV protons, 150 cm of dry air):
 * theta0 ~ 6.5 mrad over the slab, and a pencil beam ends up with a lateral
 * sigma of theta0 * t / sqrt(3) ~ 0.56 cm.  The window is wide enough to
 * survive material-table refinements but far tighter than the 0.088 cm the
 * unscattered pencil produced. */
#define SIGMA_MIN_CM 0.45
#define SIGMA_MAX_CM 0.70

/* Agreement required between the three runs (fraction of sigma).  Each run is
 * 4000 independent primaries, so the statistical spread on sigma is ~1%. */
#define SIGMA_REL_TOL 0.06

static int tmp_counter = 0;

static char const geo_one_zone[] = "    0    0  issue #325: 2 m of air as a single zone\n"
                                   "  RCC    1       0.0       0.0       0.0       0.0       0.0     200.0\n"
                                   "                 5.0\n"
                                   "  RCC    2       0.0       0.0     -10.0       0.0       0.0     215.0\n"
                                   "                 6.0\n"
                                   "  RCC    3       0.0       0.0     -11.0       0.0       0.0     217.0\n"
                                   "                 7.0\n"
                                   "  END\n"
                                   "  001          +1\n"
                                   "  002          +2     -1\n"
                                   "  003          +3     -2\n"
                                   "  END\n"
                                   "  ASSIGNMAT Air 001\n"
                                   "  ASSIGNMAT vacuum 002\n"
                                   "  ASSIGNMAT blackhole 003\n";

/* The same 2 m of air, cut into five 40 cm zones of the same material. */
static char const geo_five_zones[] = "    0    0  issue #325: the same 2 m of air split into five zones\n"
                                     "  RCC    1       0.0       0.0       0.0       0.0       0.0      40.0\n"
                                     "                 5.0\n"
                                     "  RCC    2       0.0       0.0      40.0       0.0       0.0      40.0\n"
                                     "                 5.0\n"
                                     "  RCC    3       0.0       0.0      80.0       0.0       0.0      40.0\n"
                                     "                 5.0\n"
                                     "  RCC    4       0.0       0.0     120.0       0.0       0.0      40.0\n"
                                     "                 5.0\n"
                                     "  RCC    5       0.0       0.0     160.0       0.0       0.0      40.0\n"
                                     "                 5.0\n"
                                     "  RCC    6       0.0       0.0     -10.0       0.0       0.0     215.0\n"
                                     "                 6.0\n"
                                     "  RCC    7       0.0       0.0     -11.0       0.0       0.0     217.0\n"
                                     "                 7.0\n"
                                     "  END\n"
                                     "  001          +1\n"
                                     "  002          +2\n"
                                     "  003          +3\n"
                                     "  004          +4\n"
                                     "  005          +5\n"
                                     "  006          +6     -1     -2     -3     -4     -5\n"
                                     "  007          +7     -6\n"
                                     "  END\n"
                                     "  ASSIGNMAT Air 001\n"
                                     "  ASSIGNMAT Air 002\n"
                                     "  ASSIGNMAT Air 003\n"
                                     "  ASSIGNMAT Air 004\n"
                                     "  ASSIGNMAT Air 005\n"
                                     "  ASSIGNMAT vacuum 006\n"
                                     "  ASSIGNMAT blackhole 007\n";

static char const mat_air[] = "Material Air\n"
                              "    ICRU 104\n";

static void write_temp_file(char *path_out, size_t cap, char const *content) {
    FILE *fp;

    snprintf(path_out, cap, "osh_mcs_interior_scorer_%d.tmp", tmp_counter++);
    fp = fopen(path_out, "w");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fputs(content, fp) >= 0);
    fclose(fp);
}

/*
 * Lateral sigma of the radial fluence profile written by the Cyl scorer.
 *
 * The TEXT cyl file carries one "z r_centre fluence" data line per radial bin.
 * Fluence is per unit area, so the particle count in a bin is fluence times the
 * annulus area; for a 2D Gaussian <r^2> = 2 sigma^2.
 */
static double profile_sigma_cm(char const *path) {
    FILE *fp;
    char line[4096];
    char *p;        /* cursor into `line`, past leading whitespace */
    double bin_w;   /* radial bin width [cm]                       */
    double r_inner; /* inner edge of the current annulus [cm]      */
    double r_c;     /* bin centre read back from the file [cm]     */
    double fluence; /* per-area value of the bin                   */
    double area;    /* annulus area [cm^2]                         */
    double num;     /* sum of count * r_c^2                        */
    double den;     /* sum of count                                */

    bin_w = PLANE_R_MAX_CM / (double) PLANE_NR;
    num = 0.0;
    den = 0.0;
    r_inner = 0.0;

    fp = fopen(path, "r");
    ASSERT_TRUE(fp != NULL);
    while (fgets(line, sizeof(line), fp)) {
        p = line;
        while (*p == ' ' || *p == '\t') {
            ++p;
        }
        if (*p == '#' || *p == '\n' || *p == '\r' || *p == '\0') {
            continue;
        }
        if (sscanf(p, "%*f %lf %lf", &r_c, &fluence) != 2) {
            continue;
        }
        area = OSH_M_PI * ((r_inner + bin_w) * (r_inner + bin_w) - r_inner * r_inner);
        r_inner += bin_w;
        num += fluence * area * r_c * r_c;
        den += fluence * area;
    }
    fclose(fp);

    ASSERT_TRUE(den > 0.0);
    return sqrt(0.5 * num / den);
}

/* Run one deck end to end and return the scored lateral sigma [cm]. */
static double run_case(char const *geo_text, char const *deltae) {
    char geo_path[512];
    char mat_path[512];
    char beam_path[512];
    char scoring_path[512];
    char beam_text[1024];
    char scoring_text[1024];
    char const *out_path = "osh_mcs_interior_scorer_out.dat";
    struct osh_geometry_workspace *geo = NULL;
    struct osh_beam_workspace *beam = NULL;
    struct osh_material_workspace *mat = NULL;
    struct osh_scoring_workspace *scoring = NULL;
    struct osh_simulation *sim = NULL;
    double sigma;

    snprintf(beam_text,
             sizeof(beam_text),
             "RNDSEED   12345\n"
             "PRIMARY   1 1\n"
             "TMAX0     70.2 0.0\n"
             "BEAMSIGMA 0.0 0.0\n" /* true pencil: no source width  */
             "BEAMDIV   0.0 0.0\n" /* and no source divergence      */
             "BEAMPOS   0.0 0.0 -5.0\n"
             "NSTAT     4000 -1\n"
             "DELTAE    %s\n"
             "DEMIN     0.025\n"
             "STRAGG    0\n"  /* isolate MCS: no energy fluctuation */
             "MSCAT     1\n"  /* Gaussian: exactly substep-additive */
             "NUCRE     0\n", /* no hadronic deflection channel     */
             deltae);

    snprintf(scoring_text,
             sizeof(scoring_text),
             "Geometry Cyl\n"
             "  Name Plane\n"
             "  R 0.0 %g %d\n"
             "  Z %g %g 1\n"
             "Output\n"
             "  Filename %s\n"
             "  Fileformat TEXT\n"
             "  Geo Plane\n"
             "  Quantity Fluence\n",
             PLANE_R_MAX_CM,
             PLANE_NR,
             PLANE_Z_CM - 0.1,
             PLANE_Z_CM + 0.1,
             out_path);

    write_temp_file(geo_path, sizeof(geo_path), geo_text);
    write_temp_file(mat_path, sizeof(mat_path), mat_air);
    write_temp_file(beam_path, sizeof(beam_path), beam_text);
    write_temp_file(scoring_path, sizeof(scoring_path), scoring_text);

    ASSERT_TRUE(osh_geometry_setup_from_path(geo_path, NULL, &geo) == OSH_OK);
    ASSERT_TRUE(osh_beam_setup_from_path(beam_path, NULL, &beam) == OSH_OK);
    ASSERT_TRUE(osh_material_setup_from_path(mat_path, NULL, &mat) == OSH_OK);
    ASSERT_TRUE(osh_scoring_setup_from_path(scoring_path, NULL, &scoring) == OSH_OK);

    ASSERT_TRUE(osh_simulation_create(beam, geo, mat, scoring, NULL, &sim) == OSH_OK);
    ASSERT_TRUE(osh_simulation_run(sim) == OSH_OK);
    ASSERT_TRUE(osh_simulation_save(sim) == OSH_OK);

    sigma = profile_sigma_cm(out_path);

    ASSERT_TRUE(osh_simulation_free(sim) == OSH_OK);
    osh_geometry_workspace_free(geo);
    osh_beam_workspace_free(beam);
    osh_material_workspace_free(mat);
    osh_scoring_workspace_free(scoring);
    remove(out_path);
    remove(geo_path);
    remove(mat_path);
    remove(beam_path);
    remove(scoring_path);

    return sigma;
}

static void test_interior_scorer_sees_scattering(void) {
    double one_coarse;  /* one zone, DELTAE 0.05  — the reported failure */
    double one_fine;    /* one zone, DELTAE 0.005 — fine-step reference  */
    double five_coarse; /* five zones, DELTAE 0.05 — the workaround deck */

    one_coarse = run_case(geo_one_zone, "0.05");
    one_fine = run_case(geo_one_zone, "0.005");
    five_coarse = run_case(geo_five_zones, "0.05");

    /* The bug's signature: the scorer sees an unscattered pencil, which with
     * these bins reads back as sigma ~ 0.09 cm. */
    ASSERT_TRUE(one_coarse > SIGMA_MIN_CM && one_coarse < SIGMA_MAX_CM);

    /* The scored width must not depend on the step-size control ... */
    ASSERT_TRUE(fabs(one_coarse - one_fine) < SIGMA_REL_TOL * one_fine);

    /* ... nor on whether the user subdivided the air into real zones. */
    ASSERT_TRUE(fabs(one_coarse - five_coarse) < SIGMA_REL_TOL * five_coarse);
}

int main(void) {
    test_interior_scorer_sees_scattering();
    printf("test_osh_mcs_interior_scorer: all tests passed\n");
    return 0;
}
