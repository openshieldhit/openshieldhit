#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "particle/osh_isotope_db_generated.h"
#include "particle/osh_particle.h"
#include "particle/osh_particle_const.h"
#include "particle/osh_particle_pdg.h"

#define Z_OSMIUM 76
#define EPS 1e-10

#define MASS_PROTON 938.2720882
#define MASS_NEUTRON 939.5654205
#define MASS_ELECTRON 0.51099894929
#define MASS_ANTIPROTON 938.2720882
#define MASS_DEUTERON 1875.61294257
#define MASS_HE4 3727.379378

/* Carbon-12 ion: not in particle_db, exercises the isotope_db fallback path */
#define PDG_CARBON_12 1000060120

void test_isotope_db(void) {

    struct isotope iso;

    iso = osh_isotope_db[osh_isotopes_idx_default[Z_OSMIUM]];
    assert(iso.z == 76);
    assert(iso.a == 192);
    assert((iso.symb[0] == 'O') && (iso.symb[1] == 's'));
    assert(fabs(iso.amass - 191.961477000000) < EPS);
    assert(fabs(iso.abund - 0.4078) < EPS);
}

void test_particles(void) {

    /* get a protons */
    struct particle p;

    int res;

    res = osh_particle_from_pdg(&p, OSH_PART_PDG_PROTON);
    assert(res == 1);

    assert(p.pdg == OSH_PART_PDG_PROTON);
    assert(p.charge == 1);
    assert(p.mass == MASS_PROTON);
    assert(p.is_nucleus == 0);
}

void test_particle_name_from_pdg(void) {
    char buf[64];
    int res;

    /* named particles resolved via particle_db */
    res = osh_particle_name_from_pdg(OSH_PART_PDG_PROTON, buf, sizeof(buf));
    assert(res == 1);
    assert(strcmp(buf, "proton") == 0);

    res = osh_particle_name_from_pdg(OSH_PART_PDG_DEUTERON, buf, sizeof(buf));
    assert(res == 1);
    assert(strcmp(buf, "deuteron") == 0);

    res = osh_particle_name_from_pdg(OSH_PART_PDG_HE4, buf, sizeof(buf));
    assert(res == 1);
    assert(strcmp(buf, "He-4") == 0);

    /* ion not in particle_db: falls through to isotope_db -> "C-12" */
    res = osh_particle_name_from_pdg(PDG_CARBON_12, buf, sizeof(buf));
    assert(res == 1);
    assert(strcmp(buf, "C-12") == 0);

    /* unknown PDG: returns 0 and nulls the buffer */
    buf[0] = 'x';
    res = osh_particle_name_from_pdg(999, buf, sizeof(buf));
    assert(res == 0);
    assert(buf[0] == '\0');

    /* buf_size == 0: returns 0 without touching buffer */
    res = osh_particle_name_from_pdg(OSH_PART_PDG_PROTON, buf, 0);
    assert(res == 0);
}

void test_particle_symbol_from_pdg(void) {
    char buf[64];
    int res;

    /* named particles resolved via particle_db */
    res = osh_particle_symbol_from_pdg(OSH_PART_PDG_PROTON, buf, sizeof(buf));
    assert(res == 1);
    assert(strcmp(buf, "p") == 0);

    res = osh_particle_symbol_from_pdg(OSH_PART_PDG_DEUTERON, buf, sizeof(buf));
    assert(res == 1);
    assert(strcmp(buf, "d") == 0);

    res = osh_particle_symbol_from_pdg(OSH_PART_PDG_HE4, buf, sizeof(buf));
    assert(res == 1);
    assert(strcmp(buf, "He4") == 0);

    /* ion not in particle_db: falls through to isotope_db -> "C-12" */
    res = osh_particle_symbol_from_pdg(PDG_CARBON_12, buf, sizeof(buf));
    assert(res == 1);
    assert(strcmp(buf, "C-12") == 0);

    /* unknown PDG: returns 0 and nulls the buffer */
    buf[0] = 'x';
    res = osh_particle_symbol_from_pdg(999, buf, sizeof(buf));
    assert(res == 0);
    assert(buf[0] == '\0');

    /* buf_size == 0: returns 0 without touching buffer */
    res = osh_particle_symbol_from_pdg(OSH_PART_PDG_PROTON, buf, 0);
    assert(res == 0);
}

void test_particle_from_pdg_ion(void) {
    struct particle p;
    int res;

    /* deuteron: z=1, a=2, mass from particle_db */
    res = osh_particle_from_pdg(&p, OSH_PART_PDG_DEUTERON);
    assert(res == 1);
    assert(p.pdg == OSH_PART_PDG_DEUTERON);
    assert(p.is_nucleus == 1);
    assert(p.z == 1);
    assert(p.a == 2);
    assert(p.charge == 1);
    assert(fabs(p.mass - MASS_DEUTERON) < EPS);

    /* He-4: z=2, a=4, mass from particle_db */
    res = osh_particle_from_pdg(&p, OSH_PART_PDG_HE4);
    assert(res == 1);
    assert(p.is_nucleus == 1);
    assert(p.z == 2);
    assert(p.a == 4);
    assert(p.charge == 2);
    assert(fabs(p.mass - MASS_HE4) < EPS);

    /* Carbon-12: ion path via isotope_db, nuclear mass (approximate) */
    res = osh_particle_from_pdg(&p, PDG_CARBON_12);
    assert(res == 1);
    assert(p.is_nucleus == 1);
    assert(p.z == 6);
    assert(p.a == 12);
    assert(p.charge == 6);
    assert(p.mass > 0.0);
}

void test_particle_pdg_from_name(void) {
    int pdg;
    int res;

    res = osh_particle_pdg_from_name("proton", &pdg);
    assert(res == 1);
    assert(pdg == OSH_PART_PDG_PROTON);

    res = osh_particle_pdg_from_name("P", &pdg);
    assert(res == 1);
    assert(pdg == OSH_PART_PDG_PROTON);

    res = osh_particle_pdg_from_name("antiproton", &pdg);
    assert(res == 1);
    assert(pdg == OSH_PART_PDG_APROTON);

    res = osh_particle_pdg_from_name("alpha", &pdg);
    assert(res == 1);
    assert(pdg == OSH_PART_PDG_HE4);

    res = osh_particle_pdg_from_name("no_such_particle", &pdg);
    assert(res == 0);
}

void test_particle_from_name(void) {
    struct particle p;
    int res;

    res = osh_particle_from_name(&p, "He4");
    assert(res == 1);
    assert(p.pdg == OSH_PART_PDG_HE4);
    assert(p.is_nucleus == 1);
    assert(p.z == 2);
    assert(p.a == 4);
}

int main(void) {
    printf("Running osh_particle tests...\n");

    test_isotope_db();
    test_particles();
    test_particle_name_from_pdg();
    test_particle_symbol_from_pdg();
    test_particle_from_pdg_ion();
    test_particle_pdg_from_name();
    test_particle_from_name();

    printf("All tests passed.\n");
    return 0;
}
