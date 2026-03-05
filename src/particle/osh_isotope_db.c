#include "particle/osh_isotope_db.h"

#include "particle/osh_isotope_db_generated.h"
#include "particle/osh_particle_pdg.h"

/*
 * @brief Find isotope in database by atomic number and mass number.
 * @params iso pointer to isotope struct to fill with found isotope data,
 * @params z atomic number
 * @params a nucleon number
 *
 * @returns 1 if found, 0 if not found.
 */
int osh_isotope_from_za(struct isotope *iso, unsigned int z, unsigned a) {

    int i;
    int i_max;

    if (z >= OSH_ISOTOPE_DB_NELEM) {
        return 0;
    }

    i_max = osh_isotopes_idx[z] + osh_isotopes_len[z];

    for (i = osh_isotopes_idx[z]; i < i_max; ++i) {
        if (osh_isotope_db[i].a == a) {
            *iso = osh_isotope_db[i];
            return 1;
        }
    }
    return 0;
}

int osh_isotope_from_pdg(struct isotope *iso, int pdg) {

    unsigned int z, a;

    if (pdg < OSH_PART_PDG_PROTON_ION) {
        return 0;
    }

    z = (pdg - OSH_PART_PDG_HIBASE) / 10000;
    a = (pdg - OSH_PART_PDG_HIBASE) % 10000 / 10;

    return osh_isotope_from_za(iso, z, a);
}