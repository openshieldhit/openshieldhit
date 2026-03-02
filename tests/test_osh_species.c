#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "species/osh_isotope_db_generated.h"

#define Z_OSMIUM 76
#define EPS 1e-10

void test_isotope_db() {

    struct isotope iso;

    iso = osh_isotope_db[osh_isotopes_idx_default[Z_OSMIUM]];
    assert(iso.z == 76);
    assert(iso.a == 192);
    assert((iso.symb[0] == 'O') && (iso.symb[1] == 's'));
    assert(fabs(iso.amass - 191.961477000000) < EPS);
    assert(fabs(iso.abund - 0.4078) < EPS);
}


int main(void) {
    printf("Running osh_species tests...\n");

    test_isotope_db();

    printf("All tests passed.\n");
    return 0;
}
