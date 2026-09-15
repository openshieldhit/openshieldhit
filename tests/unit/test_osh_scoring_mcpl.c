#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_app_osh.h"
#include "common/osh_step.h"
#include "mcpl.h"
#include "openshieldhit/const.h"
#include "openshieldhit/scoring.h"
#include "openshieldhit/status.h"
#include "particle/osh_particle.h"
#include "scoring/runtime/osh_scoring_accumulator.h"
#include "scoring/runtime/osh_scoring_compile.h"
#include "scoring/runtime/osh_scoring_defs.h"
#include "scoring/runtime/osh_scoring_mcpl_record.h"
#include "scoring/runtime/osh_scoring_step.h"
#include "scoring/save/osh_scoring_save.h"
#include "scoring/save/osh_scoring_save_mcpl.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

static int tmp_counter = 0;

static void assert_close(double a, double b) {
    ASSERT_TRUE(fabs(a - b) < 1.0e-9);
}

static void write_temp_file(char *path, size_t path_cap, char const *content) {
    FILE *fp;

    snprintf(path, path_cap, "osh_scoring_mcpl_test_%d.tmp", tmp_counter++);
    fp = fopen(path, "w");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fputs(content, fp) >= 0);
    ASSERT_TRUE(fclose(fp) == 0);
}

/* Parses `detect`, resolves the single Zone geometry's one zone to transport
 * zone id `zone_id` (mirroring what the app does against geo.dat — see
 * test_osh_scoring_step.c's own test_score_zone_energy_fluence_dose()), and
 * compiles. Caller owns the workspace and runtime this produces, and must
 * free both with osh_scoring_workspace_free() / osh_scoring_runtime_free(). */
static void setup_one_zone_mcpl(char const *detect,
                                size_t zone_id,
                                struct osh_scoring_workspace **ws_out,
                                struct osh_scoring_runtime *rt) {
    char path[512];
    enum osh_status rc;

    write_temp_file(path, sizeof(path), detect);
    rc = osh_scoring_setup_from_path(path, NULL, ws_out);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(*ws_out != NULL);
    ASSERT_TRUE((*ws_out)->ngeometries == 1u);
    ASSERT_TRUE((*ws_out)->geometries[0].nzone_indices == 1u);

    (*ws_out)->geometries[0].zone_indices = (size_t *) calloc(1u, sizeof(size_t));
    ASSERT_TRUE((*ws_out)->geometries[0].zone_indices != NULL);
    (*ws_out)->geometries[0].zone_indices[0] = zone_id;

    memset(rt, 0, sizeof(*rt));
    rc = osh_scoring_compile(*ws_out, NULL, rt);
    ASSERT_TRUE(rc == OSH_OK);

    remove(path);
}

static char const *k_detect = "Geometry Zone\n"
                              "    Name UpstreamPlane\n"
                              "    Zone Plane1\n"
                              "\n"
                              "Output\n"
                              "    Filename osh_scoring_mcpl_test.mcpl\n"
                              "    FileFormat MCPL\n"
                              "    Geo UpstreamPlane\n"
                              "    Quantity MCPL\n"
                              "    MaxRecords 2\n";

/* Fills a step crossing zone 5 exiting at (1,2,3) cm with direction +Z,
 * ekin_exit MeV, weight wt, generation gen. */
static void fill_step(struct step *st, double ekin_exit, double wt, uint8_t gen) {
    memset(st, 0, sizeof(*st));
    st->p[0] = 1.0;
    st->p[1] = 2.0;
    st->p[2] = 0.0;
    st->p[3] = ekin_exit + 5.0;
    st->q[0] = 1.0;
    st->q[1] = 2.0;
    st->q[2] = 3.0;
    st->q[3] = ekin_exit;
    st->v[2] = 1.0;
    st->w[2] = 1.0;
    st->ds = 3.0;
    st->rho = 1.0;
    st->wt = wt;
    st->medium = 0;
    st->zone = 5;
    st->gen = gen;
}

/* Books two protons into the MCPL page's append buffer, then confirms a third
 * is rejected once MaxRecords (2, from k_detect) is exhausted, and that the
 * booked records carry the step's exit state, weight, and generation. */
static void test_score_mcpl_books_records_and_enforces_max_records(void) {
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    struct particle part;
    struct step st;
    struct osh_scoring_page_runtime *page;
    enum osh_status rc;

    setup_one_zone_mcpl(k_detect, 5u, &ws, &rt);

    ASSERT_TRUE(rt.npages == 1u);
    page = &rt.pages[0];
    ASSERT_TRUE(page->score_kind == OSH_SCORING_SCORE_MCPL);
    ASSERT_TRUE(page->acc.mcpl_capacity == 2u);
    ASSERT_TRUE(page->acc.mcpl_records != NULL);
    ASSERT_TRUE(page->acc.mcpl_count != NULL);
    ASSERT_TRUE(*page->acc.mcpl_count == 0u);

    memset(&part, 0, sizeof(part));
    part.pdg = 2212; /* proton */
    part.mass = 938.27208816;
    part.charge = 1;
    part.z = 1u;
    part.a = 1u;

    fill_step(&st, 95.0, 0.5, 0u);
    rc = osh_scoring_score_step(
        &rt, osh_scoring_runtime_master_accumulators(&rt), osh_scoring_runtime_master_scratch(&rt), &part, &st);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(*page->acc.mcpl_count == 1u);

    fill_step(&st, 60.0, 0.25, 1u);
    rc = osh_scoring_score_step(
        &rt, osh_scoring_runtime_master_accumulators(&rt), osh_scoring_runtime_master_scratch(&rt), &part, &st);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(*page->acc.mcpl_count == 2u);

    /* MaxRecords 2 is now exhausted: a third crossing must fail loudly rather
     * than silently drop the particle or grow the buffer. */
    fill_step(&st, 10.0, 1.0, 0u);
    rc = osh_scoring_score_step(
        &rt, osh_scoring_runtime_master_accumulators(&rt), osh_scoring_runtime_master_scratch(&rt), &part, &st);
    ASSERT_TRUE(rc == OSH_ESTATE);
    ASSERT_TRUE(*page->acc.mcpl_count == 2u); /* unchanged: no partial/garbage record was added */

    assert_close(page->acc.mcpl_records[0].position[0], 1.0);
    assert_close(page->acc.mcpl_records[0].position[1], 2.0);
    assert_close(page->acc.mcpl_records[0].position[2], 3.0);
    assert_close(page->acc.mcpl_records[0].direction[2], 1.0);
    assert_close(page->acc.mcpl_records[0].ekin, 95.0);
    assert_close(page->acc.mcpl_records[0].weight, 0.5);
    ASSERT_TRUE(page->acc.mcpl_records[0].pdgcode == 2212);
    ASSERT_TRUE(page->acc.mcpl_records[0].gen == 0u);

    assert_close(page->acc.mcpl_records[1].ekin, 60.0);
    assert_close(page->acc.mcpl_records[1].weight, 0.25);
    ASSERT_TRUE(page->acc.mcpl_records[1].gen == 1u);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
}

/* End-to-end: compile, score two records, save to a real .mcpl file, then
 * read it back with the vendored MCPL reader (the same core library a
 * consumer such as Geant4's mcpl_open_file()/mcpl_read() would use) and
 * confirm every field round-trips, including the "nstat" header stat and the
 * per-particle userflags generation convention. */
static void test_save_mcpl_output_round_trips_through_mcpl_reader(void) {
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    struct particle part;
    struct step st;
    enum osh_status rc;
    mcpl_file_t f;
    mcpl_particle_t const *p;
    double nstat_stat;

    setup_one_zone_mcpl(k_detect, 5u, &ws, &rt);

    memset(&part, 0, sizeof(part));
    part.pdg = -11; /* positron, to also exercise a non-proton PDG code */
    part.mass = 0.510998950;
    part.charge = 1;

    fill_step(&st, 42.5, 0.75, 2u);
    rc = osh_scoring_score_step(
        &rt, osh_scoring_runtime_master_accumulators(&rt), osh_scoring_runtime_master_scratch(&rt), &part, &st);
    ASSERT_TRUE(rc == OSH_OK);

    remove("osh_scoring_mcpl_test.mcpl");
    rc = osh_scoring_save(ws, &rt, 1000ull);
    ASSERT_TRUE(rc == OSH_OK);

    f = mcpl_open_file("osh_scoring_mcpl_test.mcpl");
    ASSERT_TRUE(mcpl_hdr_nparticles(f) == 1u);
    ASSERT_TRUE(mcpl_hdr_has_userflags(f));
    ASSERT_TRUE(mcpl_hdr_universal_pdgcode(f) == 0); /* per-particle pdgcode, not a universal one */

    nstat_stat = mcpl_hdr_stat_sum(f, "nstat");
    assert_close(nstat_stat, 1000.0);

    p = mcpl_read(f);
    ASSERT_TRUE(p != NULL);
    assert_close(p->position[0], 1.0);
    assert_close(p->position[1], 2.0);
    assert_close(p->position[2], 3.0);
    assert_close(p->direction[2], 1.0);
    assert_close(p->ekin, 42.5);
    assert_close(p->weight, 0.75);
    ASSERT_TRUE(p->pdgcode == -11);
    ASSERT_TRUE(p->userflags == 2u); /* generation, per the header comment osh_scoring_save_mcpl_output() writes */

    p = mcpl_read(f);
    ASSERT_TRUE(p == NULL); /* exactly one record was booked */

    mcpl_close_file(f);
    remove("osh_scoring_mcpl_test.mcpl");

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
}

/* osh_scoring_compile() must reject every malformed pairing up front rather
 * than let it reach the save layer: Quantity MCPL requires both FileFormat
 * MCPL and Geometry Zone, MaxRecords must be set, and FileFormat MCPL must
 * not be paired with any other quantity. */
static void test_compile_rejects_malformed_mcpl_output(void) {
    char path[512];
    struct osh_scoring_workspace *ws;
    struct osh_scoring_runtime rt;
    enum osh_status rc;
    size_t i;
    char const *bad[] = {
        /* Quantity MCPL without FileFormat MCPL. */
        "Geometry Zone\n    Name Z\n    Zone Plane1\n\nOutput\n    Filename a\n    Geo Z\n    Quantity MCPL\n    "
        "MaxRecords 10\n",
        /* FileFormat MCPL without Quantity MCPL. */
        "Geometry Zone\n    Name Z\n    Zone Plane1\n\nOutput\n    Filename a\n    FileFormat MCPL\n    Geo Z\n    "
        "Quantity Energy\n",
        /* Quantity MCPL missing MaxRecords. */
        "Geometry Zone\n    Name Z\n    Zone Plane1\n\nOutput\n    Filename a\n    FileFormat MCPL\n    Geo Z\n    "
        "Quantity MCPL\n",
        /* Quantity MCPL on a Mesh geometry, not Zone. */
        "Geometry Mesh\n    Name M\n    X -1 1 1\n    Y -1 1 1\n    Z -1 1 1\n\nOutput\n    Filename a\n    "
        "FileFormat MCPL\n    Geo M\n    Quantity MCPL\n    MaxRecords 10\n",
        /* An MCPL output carrying a second, non-MCPL page.  MCPL is page 0 so
           this trips the "exactly one page" rule rather than the FileFormat/
           Quantity pairing rule the second case above covers. */
        "Geometry Zone\n    Name Z\n    Zone Plane1\n\nOutput\n    Filename a\n    FileFormat MCPL\n    Geo Z\n    "
        "Quantity MCPL\n    MaxRecords 10\n    Quantity Energy\n",
        /* Quantity MCPL with a differential axis: a phase-space stream has no
           binning for Diff1 to apply to. */
        "Geometry Zone\n    Name Z\n    Zone Plane1\n\nOutput\n    Filename a\n    FileFormat MCPL\n    Geo Z\n    "
        "Quantity MCPL\n    MaxRecords 10\n    Diff1 0.1 200.0 100 LOG\n",
    };

    for (i = 0; i < sizeof(bad) / sizeof(bad[0]); ++i) {
        ws = NULL;
        write_temp_file(path, sizeof(path), bad[i]);
        rc = osh_scoring_setup_from_path(path, NULL, &ws);
        ASSERT_TRUE(rc == OSH_OK);
        ASSERT_TRUE(ws != NULL);
        if (ws->ngeometries > 0u && ws->geometries[0].nzone_indices > 0u) {
            ws->geometries[0].zone_indices = (size_t *) calloc(ws->geometries[0].nzone_indices, sizeof(size_t));
            ASSERT_TRUE(ws->geometries[0].zone_indices != NULL);
        }

        memset(&rt, 0, sizeof(rt));
        rc = osh_scoring_compile(ws, NULL, &rt);
        ASSERT_TRUE(rc != OSH_OK);

        osh_scoring_runtime_free(&rt);
        osh_scoring_workspace_free(ws);
        remove(path);
    }
}

/* The MaxRecords card's own parse-time diagnostics, which sit in the app's
 * detect.dat parser rather than in osh_scoring_compile(): each of these must
 * fail the parse outright, before any scoring runtime exists. */
static void test_parse_rejects_malformed_max_records(void) {
    char path[512];
    struct osh_scoring_workspace *ws;
    enum osh_status rc;
    size_t i;
    char const *bad[] = {
        /* MaxRecords before any Quantity line, so there is no page to attach to. */
        "Geometry Zone\n    Name Z\n    Zone Plane1\n\nOutput\n    Filename a\n    FileFormat MCPL\n    Geo Z\n    "
        "MaxRecords 10\n    Quantity MCPL\n",
        /* MaxRecords with no record count. */
        "Geometry Zone\n    Name Z\n    Zone Plane1\n\nOutput\n    Filename a\n    FileFormat MCPL\n    Geo Z\n    "
        "Quantity MCPL\n    MaxRecords\n",
        /* MaxRecords 0 -- a buffer that can never hold a record. */
        "Geometry Zone\n    Name Z\n    Zone Plane1\n\nOutput\n    Filename a\n    FileFormat MCPL\n    Geo Z\n    "
        "Quantity MCPL\n    MaxRecords 0\n",
    };

    for (i = 0; i < sizeof(bad) / sizeof(bad[0]); ++i) {
        ws = NULL;
        write_temp_file(path, sizeof(path), bad[i]);
        rc = osh_scoring_setup_from_path(path, NULL, &ws);
        ASSERT_TRUE(rc != OSH_OK);
        osh_scoring_workspace_free(ws);
        remove(path);
    }
}

/* osh_scoring_save_mcpl_output() re-checks the shape osh_scoring_compile()
 * already guarantees rather than trusting its caller, so those guards are
 * unreachable through the normal path.  Reach them here by handing it bad
 * arguments and by poking a compiled runtime, so the contract in the header
 * is actually exercised instead of merely asserted. */
static void test_save_mcpl_output_rejects_bad_arguments(void) {
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    size_t saved_count;
    enum osh_status rc;

    ASSERT_TRUE(osh_scoring_save_mcpl_output(NULL, NULL, 1ull, 0u) == OSH_EINVAL);

    setup_one_zone_mcpl(k_detect, 5u, &ws, &rt);
    ASSERT_TRUE(rt.noutputs == 1u);
    ASSERT_TRUE(rt.npages == 1u);

    /* output_idx past the end. */
    rc = osh_scoring_save_mcpl_output(ws, &rt, 1ull, rt.noutputs);
    ASSERT_TRUE(rc == OSH_EINVAL);

    /* nstat == 0 would make the per-primary normalisation meaningless. */
    rc = osh_scoring_save_mcpl_output(ws, &rt, 0ull, 0u);
    ASSERT_TRUE(rc == OSH_EINVAL);

    /* A multi-page output: rejected as unsupported rather than writing page 0
     * and silently dropping the rest. */
    rt.outputs[0].npages = 2u;
    rc = osh_scoring_save_mcpl_output(ws, &rt, 1ull, 0u);
    ASSERT_TRUE(rc == OSH_ENOTSUP);
    rt.outputs[0].npages = 1u;

    /* A page that is not an MCPL page at all. */
    rt.pages[0].score_kind = OSH_SCORING_SCORE_ENERGY;
    rc = osh_scoring_save_mcpl_output(ws, &rt, 1ull, 0u);
    ASSERT_TRUE(rc == OSH_ESTATE);
    rt.pages[0].score_kind = OSH_SCORING_SCORE_MCPL;

    /* A count past the pre-allocated capacity: a hot-path bound-check bug, so
     * refuse rather than read off the end of the buffer. */
    saved_count = *rt.pages[0].acc.mcpl_count;
    *rt.pages[0].acc.mcpl_count = rt.pages[0].acc.mcpl_capacity + 1u;
    rc = osh_scoring_save_mcpl_output(ws, &rt, 1ull, 0u);
    ASSERT_TRUE(rc == OSH_ESTATE);
    *rt.pages[0].acc.mcpl_count = saved_count;

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
}

/* mcpl_add_particle() aborts the whole process on a direction vector that is
 * not unit length, so the save layer normalises every record on the way out.
 * A degenerate (zero) direction cannot be divided back to unit length, and
 * maps to +Z instead.  st->w is always a unit vector in practice, so book a
 * record normally and then zero its direction to reach the fallback. */
static void test_save_mcpl_output_maps_degenerate_direction_to_z(void) {
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    struct particle part;
    struct step st;
    struct osh_scoring_page_runtime *page;
    enum osh_status rc;
    mcpl_file_t f;
    mcpl_particle_t const *p;

    setup_one_zone_mcpl(k_detect, 5u, &ws, &rt);
    page = &rt.pages[0];

    memset(&part, 0, sizeof(part));
    part.pdg = 2212;
    part.mass = 938.27208816;
    part.charge = 1;
    part.z = 1u;
    part.a = 1u;

    fill_step(&st, 30.0, 1.0, 0u);
    rc = osh_scoring_score_step(
        &rt, osh_scoring_runtime_master_accumulators(&rt), osh_scoring_runtime_master_scratch(&rt), &part, &st);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(*page->acc.mcpl_count == 1u);

    page->acc.mcpl_records[0].direction[0] = 0.0;
    page->acc.mcpl_records[0].direction[1] = 0.0;
    page->acc.mcpl_records[0].direction[2] = 0.0;

    remove("osh_scoring_mcpl_test.mcpl");
    rc = osh_scoring_save(ws, &rt, 1ull);
    ASSERT_TRUE(rc == OSH_OK);

    f = mcpl_open_file("osh_scoring_mcpl_test.mcpl");
    p = mcpl_read(f);
    ASSERT_TRUE(p != NULL);
    assert_close(p->direction[0], 0.0);
    assert_close(p->direction[1], 0.0);
    assert_close(p->direction[2], 1.0);
    mcpl_close_file(f);
    remove("osh_scoring_mcpl_test.mcpl");

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
}

/* APPEND's cross-run combine rule is concatenation, which the additive folds
 * in osh_scoring_accumulator_merge() cannot express.  Merging an accumulator
 * that carries an MCPL buffer must therefore refuse rather than silently drop
 * one side's records -- today no caller does it, but a future parallel worker
 * would. */
static void test_accumulator_merge_refuses_mcpl_records(void) {
    struct osh_scoring_accumulator dst;
    struct osh_scoring_accumulator src;
    struct osh_scoring_mcpl_record rec;
    double dst_data[1];
    double src_data[1];
    size_t count = 0u;

    memset(&dst, 0, sizeof(dst));
    memset(&src, 0, sizeof(src));
    memset(&rec, 0, sizeof(rec));
    dst_data[0] = 1.0;
    src_data[0] = 2.0;
    dst.data = dst_data;
    src.data = src_data;
    dst.len = 1u;
    src.len = 1u;

    /* Baseline: without an MCPL buffer the same pair merges fine, so the
     * refusal below is about the records and not about the shape. */
    ASSERT_TRUE(osh_scoring_accumulator_merge(&dst, &src) == OSH_OK);

    src.mcpl_records = &rec;
    src.mcpl_count = &count;
    src.mcpl_capacity = 1u;
    ASSERT_TRUE(osh_scoring_accumulator_merge(&dst, &src) == OSH_ENOTSUP);

    src.mcpl_records = NULL;
    src.mcpl_count = NULL;
    src.mcpl_capacity = 0u;
    dst.mcpl_records = &rec;
    dst.mcpl_count = &count;
    dst.mcpl_capacity = 1u;
    ASSERT_TRUE(osh_scoring_accumulator_merge(&dst, &src) == OSH_ENOTSUP);
}

/* osh_scoring_estimate_memory() has to account for the MaxRecords append
 * buffer on top of the usual bins-sized accumulator, since that buffer is the
 * dominant allocation of an MCPL page and the run-control uses this estimate
 * to decide whether a setup fits. */
static void test_estimate_memory_counts_the_mcpl_append_buffer(void) {
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_mem_estimate est;
    struct osh_scoring_mem_estimate est_no_mcpl;
    char path[512];
    enum osh_status rc;
    /* k_detect's MaxRecords 2, with the Quantity swapped for a plain Energy
     * page over the same single-zone geometry: same shape, no append buffer. */
    char const *detect_energy = "Geometry Zone\n"
                                "    Name UpstreamPlane\n"
                                "    Zone Plane1\n"
                                "\n"
                                "Output\n"
                                "    Filename osh_scoring_mcpl_test.dat\n"
                                "    FileFormat TEXT\n"
                                "    Geo UpstreamPlane\n"
                                "    Quantity Energy\n";

    write_temp_file(path, sizeof(path), k_detect);
    rc = osh_scoring_setup_from_path(path, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    memset(&est, 0, sizeof(est));
    rc = osh_scoring_estimate_memory(ws, &est);
    ASSERT_TRUE(rc == OSH_OK);
    osh_scoring_workspace_free(ws);
    remove(path);

    ws = NULL;
    write_temp_file(path, sizeof(path), detect_energy);
    rc = osh_scoring_setup_from_path(path, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    memset(&est_no_mcpl, 0, sizeof(est_no_mcpl));
    rc = osh_scoring_estimate_memory(ws, &est_no_mcpl);
    ASSERT_TRUE(rc == OSH_OK);
    osh_scoring_workspace_free(ws);
    remove(path);

    /* The MCPL page must cost exactly MaxRecords (2) extra records. */
    ASSERT_TRUE(est.accum_bytes == est_no_mcpl.accum_bytes + 2u * (uint64_t) sizeof(struct osh_scoring_mcpl_record));
}

/* The MCPL hot path re-checks invariants osh_scoring_compile() establishes, so
 * a future change that breaks one fails loudly instead of writing through a
 * null pointer or past the end of the append buffer.  Reach those guards by
 * corrupting a compiled runtime, since nothing else can produce this state. */
static void test_score_mcpl_guards_a_corrupt_accumulator(void) {
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    struct osh_scoring_accumulator *accs;
    struct osh_scoring_mcpl_record *saved_records;
    size_t saved_stride;
    struct particle part;
    struct step st;
    enum osh_status rc;

    setup_one_zone_mcpl(k_detect, 5u, &ws, &rt);
    accs = osh_scoring_runtime_master_accumulators(&rt);

    memset(&part, 0, sizeof(part));
    part.pdg = 2212;
    part.mass = 938.27208816;
    part.charge = 1;
    part.z = 1u;
    part.a = 1u;
    fill_step(&st, 30.0, 1.0, 0u);

    /* A crossing index past the page's stride. */
    saved_stride = rt.pages[0].diff_stride;
    rt.pages[0].diff_stride = 0u;
    rc = osh_scoring_score_step(&rt, accs, osh_scoring_runtime_master_scratch(&rt), &part, &st);
    ASSERT_TRUE(rc == OSH_ESTATE);
    rt.pages[0].diff_stride = saved_stride;

    /* An accumulator with no append buffer behind an MCPL page. */
    saved_records = accs[0].mcpl_records;
    accs[0].mcpl_records = NULL;
    rc = osh_scoring_score_step(&rt, accs, osh_scoring_runtime_master_scratch(&rt), &part, &st);
    ASSERT_TRUE(rc == OSH_ESTATE);
    accs[0].mcpl_records = saved_records;

    /* Booking still works once the runtime is intact again. */
    rc = osh_scoring_score_step(&rt, accs, osh_scoring_runtime_master_scratch(&rt), &part, &st);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(*rt.pages[0].acc.mcpl_count == 1u);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
}

/* A filter on a Quantity MCPL line selects which particles reach the phase
 * space, exactly as it does for any other quantity (documented in
 * docs/user/detect.dat.md).  A particle the filter rejects must book no
 * record at all rather than one with a zeroed or partial payload. */
static void test_score_mcpl_honours_page_filters(void) {
    char const *detect_filtered = "Filter\n"
                                  "    Name NeverMatches\n"
                                  "    Z = 99\n"
                                  "\n"
                                  "Geometry Zone\n"
                                  "    Name UpstreamPlane\n"
                                  "    Zone Plane1\n"
                                  "\n"
                                  "Output\n"
                                  "    Filename osh_scoring_mcpl_test.mcpl\n"
                                  "    FileFormat MCPL\n"
                                  "    Geo UpstreamPlane\n"
                                  "    Quantity MCPL NeverMatches\n"
                                  "    MaxRecords 2\n";
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    struct particle part;
    struct step st;
    enum osh_status rc;

    setup_one_zone_mcpl(detect_filtered, 5u, &ws, &rt);

    memset(&part, 0, sizeof(part));
    part.pdg = 2212;
    part.mass = 938.27208816;
    part.charge = 1;
    part.z = 1u; /* not Z = 99, so NeverMatches rejects it */
    part.a = 1u;

    fill_step(&st, 30.0, 1.0, 0u);
    rc = osh_scoring_score_step(
        &rt, osh_scoring_runtime_master_accumulators(&rt), osh_scoring_runtime_master_scratch(&rt), &part, &st);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(*rt.pages[0].acc.mcpl_count == 0u);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
}

/* The MCPL append buffer lives only on the master accumulator set: the private
 * sets that variance batching and --score-replicas deposit into are built by
 * osh_scoring_runtime_alloc_accumulator_set(), which allocates the binned
 * arrays but no record buffer, and osh_scoring_accumulator_merge() refuses to
 * concatenate records anyway.  Both combinations must therefore be refused
 * before transport starts, with a message that says why -- not discovered as a
 * bare OSH_ESTATE out of the hot path on the first crossing.
 *
 * The --score-replicas half of the pair lives in osh_simulation_run(), which
 * needs a whole simulation to exercise; this covers the detect.dat-only half
 * that osh_scoring_compile() owns, plus the helper both guards share. */
static void test_compile_rejects_mcpl_with_variance(void) {
    char path[512];
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    enum osh_status rc;
    char const *detect_variance = "Settings\n"
                                  "    Name withErr\n"
                                  "    Variance On\n"
                                  "\n"
                                  "Geometry Zone\n"
                                  "    Name UpstreamPlane\n"
                                  "    Zone Plane1\n"
                                  "\n"
                                  "Output\n"
                                  "    Filename osh_scoring_mcpl_test.dat\n"
                                  "    FileFormat TEXT\n"
                                  "    Geo UpstreamPlane\n"
                                  "    Quantity Energy withErr\n"
                                  "\n"
                                  "Output\n"
                                  "    Filename osh_scoring_mcpl_test.mcpl\n"
                                  "    FileFormat MCPL\n"
                                  "    Geo UpstreamPlane\n"
                                  "    Quantity MCPL\n"
                                  "    MaxRecords 10\n";

    /* Variance is run-wide: it is enabled on the Energy page, not the MCPL one,
       and must still be caught. */
    write_temp_file(path, sizeof(path), detect_variance);
    rc = osh_scoring_setup_from_path(path, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(ws != NULL);
    ws->geometries[0].zone_indices = (size_t *) calloc(1u, sizeof(size_t));
    ASSERT_TRUE(ws->geometries[0].zone_indices != NULL);
    ws->geometries[0].zone_indices[0] = 5u;

    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_ENOTSUP);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(path);

    /* The same detect.dat without Variance On compiles, so the rejection above
       is about the pairing and not about the two-output shape. */
    ws = NULL;
    setup_one_zone_mcpl(k_detect, 5u, &ws, &rt);
    ASSERT_TRUE(osh_scoring_runtime_has_mcpl_page(&rt) == 1);
    ASSERT_TRUE(osh_scoring_runtime_tracks_variance(&rt) == 0);
    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
}

/* The helper both guards are built on must not fire on a runtime with no MCPL
 * page at all, or every variance/replica run would be refused. */
static void test_has_mcpl_page_is_false_without_an_mcpl_output(void) {
    char path[512];
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    enum osh_status rc;
    char const *detect_energy = "Geometry Zone\n"
                                "    Name UpstreamPlane\n"
                                "    Zone Plane1\n"
                                "\n"
                                "Output\n"
                                "    Filename osh_scoring_mcpl_test.dat\n"
                                "    FileFormat TEXT\n"
                                "    Geo UpstreamPlane\n"
                                "    Quantity Energy\n";

    ASSERT_TRUE(osh_scoring_runtime_has_mcpl_page(NULL) == 0);

    write_temp_file(path, sizeof(path), detect_energy);
    rc = osh_scoring_setup_from_path(path, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    ws->geometries[0].zone_indices = (size_t *) calloc(1u, sizeof(size_t));
    ASSERT_TRUE(ws->geometries[0].zone_indices != NULL);

    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(osh_scoring_runtime_has_mcpl_page(&rt) == 0);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(path);
}

/* A MaxRecords overflow surfaces at the transport call site as a bare
 * OSH_ESTATE, which is indistinguishable there from an internal-invariant
 * violation and names neither MCPL nor the card to raise.  The failure path
 * calls osh_scoring_runtime_mcpl_full_page() to attribute it: NULL while the
 * buffer still has room, the offending page once it is full, with the output
 * filename and the limit reachable from it. */
static void test_mcpl_full_page_attributes_a_maxrecords_overflow(void) {
    struct osh_scoring_workspace *ws = NULL;
    struct osh_scoring_runtime rt;
    struct osh_scoring_page_runtime const *full;
    struct particle part;
    struct step st;
    enum osh_status rc;
    int i;

    ASSERT_TRUE(osh_scoring_runtime_mcpl_full_page(NULL) == NULL);

    setup_one_zone_mcpl(k_detect, 5u, &ws, &rt); /* MaxRecords 2 */
    ASSERT_TRUE(osh_scoring_runtime_mcpl_full_page(&rt) == NULL);

    memset(&part, 0, sizeof(part));
    part.pdg = 2212;
    part.mass = 938.27208816;
    part.charge = 1;
    part.z = 1u;
    part.a = 1u;

    for (i = 0; i < 2; ++i) {
        fill_step(&st, 95.0, 1.0, 0u);
        rc = osh_scoring_score_step(
            &rt, osh_scoring_runtime_master_accumulators(&rt), osh_scoring_runtime_master_scratch(&rt), &part, &st);
        ASSERT_TRUE(rc == OSH_OK);
    }
    /* Reported as full as soon as the last slot is taken, i.e. before the
     * crossing that actually fails: the attribution has to survive being asked
     * after the run has already unwound. */
    full = osh_scoring_runtime_mcpl_full_page(&rt);
    ASSERT_TRUE(full == &rt.pages[0]);
    ASSERT_TRUE(full->acc.mcpl_capacity == 2u);
    ASSERT_TRUE(*full->acc.mcpl_count == 2u);
    ASSERT_TRUE(full->output_idx < rt.noutputs);
    ASSERT_TRUE(strcmp(rt.outputs[full->output_idx].filename, "osh_scoring_mcpl_test.mcpl") == 0);

    fill_step(&st, 10.0, 1.0, 0u);
    rc = osh_scoring_score_step(
        &rt, osh_scoring_runtime_master_accumulators(&rt), osh_scoring_runtime_master_scratch(&rt), &part, &st);
    ASSERT_TRUE(rc == OSH_ESTATE);
    ASSERT_TRUE(osh_scoring_runtime_mcpl_full_page(&rt) == &rt.pages[0]);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
}

int main(void) {
    test_score_mcpl_books_records_and_enforces_max_records();
    test_save_mcpl_output_round_trips_through_mcpl_reader();
    test_compile_rejects_malformed_mcpl_output();
    test_parse_rejects_malformed_max_records();
    test_save_mcpl_output_rejects_bad_arguments();
    test_save_mcpl_output_maps_degenerate_direction_to_z();
    test_accumulator_merge_refuses_mcpl_records();
    test_estimate_memory_counts_the_mcpl_append_buffer();
    test_score_mcpl_guards_a_corrupt_accumulator();
    test_score_mcpl_honours_page_filters();
    test_compile_rejects_mcpl_with_variance();
    test_has_mcpl_page_is_false_without_an_mcpl_output();
    test_mcpl_full_page_attributes_a_maxrecords_overflow();
    printf("All MCPL scoring tests passed.\n");
    return 0;
}
