/*
 * App-layer coverage for the periodic-dump wiring in osh_run() (issue #193):
 * the time-cadence report, the beam.dat-sourced cadence report labels, the
 * shadow memory-budget reservation, and the over-budget refusal.  Drives
 * osh_run() end to end against the bundled 12_partial_dump water inputs
 * (geo/beam/mat/Water.txt), overriding only the detect file with a DoseGy scorer
 * so the snapshot shadow has a non-zero size (DoseGy postprocess writes data,
 * unlike Energy/Fluence), which is what makes the reservation branch fire.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_run.h"
#include "common/osh_file.h"
#include "openshieldhit/status.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

static int tmp_counter = 0;

static void write_temp_file(char *path_out, size_t cap, char const *content) {
    FILE *fp;
    snprintf(path_out, cap, "osh_rundump_%d.tmp", tmp_counter++);
    fp = fopen(path_out, "w");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fputs(content, fp) >= 0);
    ASSERT_TRUE(fclose(fp) == 0);
}

/* A single-quantity DoseGy mesh; 200 bins → ~1.6 KiB accumulator + ~1.6 KiB
 * snapshot shadow, so the reservation is measurable and a 1 KiB budget refuses. */
static char const *const DOSEGY_DETECT = "Geometry Mesh\n"
                                         "  Name M\n"
                                         "  X -5.0 5.0 1\n"
                                         "  Y -5.0 5.0 1\n"
                                         "  Z 0.0 20.0 200\n"
                                         "Output\n"
                                         "  Filename dose.dat\n"
                                         "  Fileformat TEXT\n"
                                         "  Geo M\n"
                                         "  Quantity DoseGy\n";

static void init_opts(struct osh_run_options *opt, char const *detect_path, char const *out_dir) {
    memset(opt, 0, sizeof(*opt));
    opt->workdir = OSH_PROJECT_SOURCE_DIR "/tests/cases/12_partial_dump";
    opt->detect_path = detect_path; /* override just the scorer; geo/beam/mat come from workdir */
    opt->out_dir = out_dir;
    opt->has_nstat = 1;
    opt->nstat = 40ull; /* tiny: the branches under test are all at setup time */
}

/*
 * A scheduled time-cadence dump run: exercises the "Periodic dumps: every N s"
 * report and the shadow-reservation accounting/report (the DoseGy page makes
 * shadow_bytes > 0).  A non-NULL out stream is required so the reporting
 * branches are taken.
 */
static void test_time_cadence_run_reports_and_reserves(void) {
    char detect_path[256];
    char out_dir[256];
    struct osh_run_options opt;
    FILE *out;
    FILE *err;

    write_temp_file(detect_path, sizeof(detect_path), DOSEGY_DETECT);
    snprintf(out_dir, sizeof(out_dir), "osh_rundump_out_%d", tmp_counter++);

    init_opts(&opt, detect_path, out_dir);
    opt.has_dump_every = 1;
    opt.dump_every_s = 0.05; /* time cadence: reported + reserved regardless of whether a dump fires */

    out = tmpfile();
    err = tmpfile();
    ASSERT_TRUE(out != NULL && err != NULL);
    ASSERT_TRUE(osh_run(&opt, out, err) == OSH_OK);
    fclose(out);
    fclose(err);
    remove(detect_path);
    osh_path_remove_dir(out_dir);
}

/* A valid beam.dat whose cadence comes from the cards themselves — a positive
 * NSTAT save step and a DUMPEVERY duration — with no CLI override.  Driving
 * osh_run() with this exercises the beam.dat-sourced report labels ("beam.dat
 * NSTAT step" / "beam.dat DUMPEVERY"), the else arms of the CLI-over-card
 * precedence, which the CLI-flag tests above never reach. */
static char const *const BEAM_WITH_CARDS = "RNDSEED   89736501\n"
                                           "PRIMARY   1 1\n"
                                           "TMAX0     150.0 0.0\n"
                                           "BEAMSIGMA 0.0 0.0\n"
                                           "BEAMPOS   0.0 0.0 -5.00\n"
                                           "NSTAT     200 40\n" /* save step 40 -> nsave = 40 */
                                           "DUMPEVERY 0.05\n"   /* wall-time cadence from the card */
                                           "DELTAE    0.005\n"
                                           "DEMIN     0.025\n"
                                           "STRAGG    1\n"
                                           "MSCAT     2\n"
                                           "NUCRE     0\n";

/*
 * Cadence sourced from beam.dat rather than the CLI: with no --dump-every[-primaries]
 * flag but a beam.dat carrying a positive NSTAT step and a DUMPEVERY duration, both
 * "Periodic dumps" report lines take their else arm and print the "beam.dat ..."
 * source label.  A non-NULL out stream is required so the reporting branches run.
 */
static void test_beam_dat_cadence_reports(void) {
    char detect_path[256];
    char beam_path[256];
    char out_dir[256];
    struct osh_run_options opt;
    FILE *out;
    FILE *err;

    write_temp_file(detect_path, sizeof(detect_path), DOSEGY_DETECT);
    write_temp_file(beam_path, sizeof(beam_path), BEAM_WITH_CARDS);
    snprintf(out_dir, sizeof(out_dir), "osh_rundump_out_%d", tmp_counter++);

    init_opts(&opt, detect_path, out_dir);
    opt.beam_path = beam_path; /* override beam.dat; no has_dump_every* -> cadence comes from the cards */

    out = tmpfile();
    err = tmpfile();
    ASSERT_TRUE(out != NULL && err != NULL);
    ASSERT_TRUE(osh_run(&opt, out, err) == OSH_OK);
    fclose(out);
    fclose(err);
    remove(beam_path);
    remove(detect_path);
    osh_path_remove_dir(out_dir);
}

/*
 * The budget-reservation refusal: with a scheduled dump the snapshot shadow is
 * added to the footprint, and a tiny --mem-budget must refuse the run up front
 * (OSH_ENOMEM) before any transport, exercising the "drop periodic dumps" hint.
 * A non-NULL err stream is required so the refusal message is written.
 */
static void test_scheduled_dump_over_budget_is_refused(void) {
    char detect_path[256];
    char out_dir[256];
    struct osh_run_options opt;
    FILE *err;

    write_temp_file(detect_path, sizeof(detect_path), DOSEGY_DETECT);
    snprintf(out_dir, sizeof(out_dir), "osh_rundump_out_%d", tmp_counter++);

    init_opts(&opt, detect_path, out_dir);
    opt.has_dump_every_primaries = 1;
    opt.dump_every_primaries = 40ull;
    opt.mem_budget = "1KB"; /* accumulator + reserved shadow (~3 KiB) exceeds this */

    err = tmpfile();
    ASSERT_TRUE(err != NULL);
    ASSERT_TRUE(osh_run(&opt, NULL, err) == OSH_ENOMEM);
    fclose(err);
    remove(detect_path);
    osh_path_remove_dir(out_dir);
}

int main(void) {
    test_time_cadence_run_reports_and_reserves();
    test_beam_dat_cadence_reports();
    test_scheduled_dump_over_budget_is_refused();
    return 0;
}
