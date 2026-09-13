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
#include "scoring/runtime/osh_scoring_compile.h"
#include "scoring/runtime/osh_scoring_defs.h"
#include "scoring/runtime/osh_scoring_mcpl_record.h"
#include "scoring/runtime/osh_scoring_step.h"
#include "scoring/save/osh_scoring_save.h"

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

int main(void) {
    test_score_mcpl_books_records_and_enforces_max_records();
    test_save_mcpl_output_round_trips_through_mcpl_reader();
    test_compile_rejects_malformed_mcpl_output();
    printf("All MCPL scoring tests passed.\n");
    return 0;
}
