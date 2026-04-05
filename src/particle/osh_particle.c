#include "particle/osh_particle.h"

#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include "common/osh_const.h"
#include "common/osh_logger.h"
#include "particle/osh_isotope_db.h"
#include "particle/osh_particle_const.h"
#include "particle/osh_particle_db.h"
#include "particle/osh_particle_pdg.h"

static void reset_particle(struct particle *p) {
    p->mass = 0.0;
    p->weight = 1.0;
    p->gen = 0;
    p->nprim = 0;
    p->pdg = OSH_PART_PDG_NONE;
    p->charge = 0;
    p->z = 0;
    p->a = 0;
    p->is_nucleus = 0;
}

static int _name_eq(char const *lhs, char const *rhs) {
    unsigned char cl;
    unsigned char cr;

    if (!lhs || !rhs) {
        return 0;
    }

    while (*lhs != '\0' && *rhs != '\0') {
        cl = (unsigned char) tolower((unsigned char) *lhs);
        cr = (unsigned char) tolower((unsigned char) *rhs);
        if (cl != cr) {
            return 0;
        }
        lhs++;
        rhs++;
    }
    return *lhs == '\0' && *rhs == '\0';
}

int osh_particle_from_pdg(struct particle *p, int pdg) {

    size_t i;
    struct isotope iso;

    p->gen = 0;
    p->nprim = 0;
    p->weight = 1.0;

    if (pdg == OSH_PART_PDG_PROTON_ION) { /* hypothetical hydrogen nucleus */
        pdg = OSH_PART_PDG_PROTON;
    }

    if (pdg == OSH_PART_PDG_NEUTRON_ION) { /* hypothetical neutron nucleus */
        pdg = OSH_PART_PDG_NEUTRON;
    }

    /* set the mass of the particle */
    if (!osh_particle_mass_from_pdg(pdg, &p->mass)) {
        reset_particle(p);
        p->pdg = OSH_PART_PDG_INVALID;
        return 0;
    }

    if (osh_particle_pdg_is_ion(pdg)) {
        if (!osh_isotope_from_pdg(&iso, pdg)) {
            reset_particle(p);
            p->pdg = OSH_PART_PDG_INVALID;
            return 0;
        }
        p->pdg = pdg;
        p->is_nucleus = 1;
        p->z = iso.z;
        p->a = iso.a;
        p->charge = (int16_t) iso.z; /* assume fully ionized */
        return 1;
    } else {
        p->is_nucleus = 0;
        p->z = 0;
        p->a = 0;
        for (i = 0; i < osh_particle_db_len; ++i) {
            if (osh_particle_db[i].pdg == pdg) {
                p->charge = osh_particle_db[i].charge_e;
                p->pdg = pdg;
                return 1;
            }
        }
    }
    reset_particle(p);
    p->pdg = OSH_PART_PDG_INVALID;
    return 0;
}

int osh_particle_pdg_from_name(char const *name, int *pdg_out) {
    size_t i;

    if (!name || !pdg_out) {
        return 0;
    }

    for (i = 0; i < osh_particle_db_len; ++i) {
        if (_name_eq(name, osh_particle_db[i].name) || _name_eq(name, osh_particle_db[i].symbol)) {
            *pdg_out = osh_particle_db[i].pdg;
            return 1;
        }
    }

    if (_name_eq(name, "alpha")) {
        *pdg_out = OSH_PART_PDG_HE4;
        return 1;
    }
    if (_name_eq(name, "pbar")) {
        *pdg_out = OSH_PART_PDG_APROTON;
        return 1;
    }
    if (_name_eq(name, "nbar")) {
        *pdg_out = OSH_PART_PDG_ANEUTRON;
        return 1;
    }

    return 0;
}

int osh_particle_from_name(struct particle *p, char const *name) {
    int pdg;

    if (!osh_particle_pdg_from_name(name, &pdg)) {
        return 0;
    }
    return osh_particle_from_pdg(p, pdg);
}

int osh_particle_pdg_is_ion(int pdg) {
    return pdg > OSH_PART_PDG_HIBASE;
}

int osh_particle_name_from_pdg(int pdg, char *name_buf, size_t buf_size) {
    struct isotope iso;
    size_t i;

    if (buf_size == 0) {
        return 0;
    }

    for (i = 0; i < osh_particle_db_len; ++i) {
        if (osh_particle_db[i].pdg == pdg) {
            snprintf(name_buf, buf_size, "%s", osh_particle_db[i].name);
            return 1;
        }
    }

    /* Isotope names are not stored in particle db, so look up in isotope database for now. */
    if (osh_isotope_from_pdg(&iso, pdg)) {
        snprintf(name_buf, buf_size, "%s-%u", iso.symb, iso.a);
        return 1;
    }
    name_buf[0] = '\0';
    return 0; /* not found */
}

int osh_particle_symbol_from_pdg(int pdg, char *symbol_buf, size_t buf_size) {
    struct isotope iso;
    size_t i;

    if (buf_size == 0) {
        return 0;
    }

    for (i = 0; i < osh_particle_db_len; ++i) {
        if (osh_particle_db[i].pdg == pdg) {
            snprintf(symbol_buf, buf_size, "%s", osh_particle_db[i].symbol);
            return 1;
        }
    }

    /* try looking up in isotope database */
    if (osh_isotope_from_pdg(&iso, pdg)) {
        snprintf(symbol_buf, buf_size, "%s-%u", iso.symb, iso.a);
        return 1;
    }
    symbol_buf[0] = '\0';
    return 0; /* not found */
}

int osh_particle_mass_from_pdg(int pdg, double *mass) {

    size_t i;
    struct isotope iso;

    /* check first if particle is in particle pdg list */
    for (i = 0; i < osh_particle_db_len; ++i) {
        if (osh_particle_db[i].pdg == pdg) {
            *mass = osh_particle_db[i].mass_mev;
            return 1;
        }
    }

    /* For ions with z = 1 and 2 always prefer masses from PDG list, if available
       since these are more accurate and include electron binding energy, which is important for light ions.
       If not found, check if it is an ion and look up in isotope db, with approximate nuclear mass.
    */
    if (osh_particle_pdg_is_ion(pdg)) {
        if (osh_isotope_from_pdg(&iso, pdg)) {
            /* convert from atomic mass to nuclear mass */
            *mass = iso.amass * OSH_AMU - (double) iso.z * OSH_PART_MASS_ELECTRON;
            return 1;
        }
    }

    return 0; /* not found */
}

void osh_print_particle(struct particle const *p) {
    char name_buf[64];

    if (!p) {
        return;
    }
    if (!osh_particle_name_from_pdg(p->pdg, name_buf, sizeof(name_buf))) {
        name_buf[0] = '\0';
    }

    osh_info("Particle: %s", name_buf);
    osh_info(OSH_LOG_HLINE);
    osh_info("%-18s : %i", "PDG code", p->pdg);
    osh_info("%-18s : %i", "Z", (int) p->z);
    osh_info("%-18s : %i", "A", (int) p->a);
    /* Mass is the bare nuclear mass (fully stripped ion, no electrons).
     * Derived as: atomic_mass * OSH_AMU - Z * m_electron (CODATA 2018). */
    osh_info("%-18s : %-12.5f MeV/c^2  (nuclear, CODATA 2018)", "Mass", p->mass);
    osh_info("%-18s : %-12.5f Da        (nuclear, CODATA 2018)", "Mass", p->mass / OSH_AMU);
    osh_info("%-18s : %i e", "Charge", (int) p->charge);
    osh_info("%-18s : %u", "Generation", p->gen);
    osh_info("%-18s : %f", "Stat.weight", p->weight);
}
