#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "particle/osh_isotope_db_generated.h"
#include "particle/osh_particle.h"
#include "particle/osh_particle_const.h"

#define Z_OSMIUM 76
#define EPS 1e-10

#define MASS_PROTON 938.2720882
#define MASS_NEUTRON 939.5654205
#define MASS_ELECTRON 0.51099894929
#define MASS_ANTIPROTON 938.2720882

void test_isotope_db() {

    struct isotope iso;

    iso = osh_isotope_db[osh_isotopes_idx_default[Z_OSMIUM]];
    assert(iso.z == 76);
    assert(iso.a == 192);
    assert((iso.symb[0] == 'O') && (iso.symb[1] == 's'));
    assert(fabs(iso.amass - 191.961477000000) < EPS);
    assert(fabs(iso.abund - 0.4078) < EPS);
}

void test_particles() {

    /* get a protons */
    struct particle p;

    int res = osh_particle_from_pdg(&p, OSH_PART_PDG_PROTON);
    assert(res == 1);

    assert(p.pdg == OSH_PART_PDG_PROTON);
    assert(p.charge == 1);
    assert(p.mass == MASS_PROTON);
    assert(p.is_nucleus == 0);
}

int main(void) {
    printf("Running osh_particle tests...\n");

    test_isotope_db();
    test_particles();

    printf("All tests passed.\n");
    return 0;
}
