#include "particle/osh_particle.h"

#include <stdio.h>

#include "particle/osh_isotope_db.h"
#include "particle/osh_particle_const.h"
#include "particle/osh_particle_db.h"
#include "particle/osh_particle_pdg.h"

int osh_particle_from_pdg(struct particle *p, int pdg) {

    size_t i;
    struct isotope iso;

    p->gen = 0;
    p->nprim = 0;
    p->weight = 1.0;

    if (pdg == 1000010010)
        pdg = OSH_PART_PDG_PROTON;

    if (pdg == 1000000010) /* hypothetical neutron nucleus */
        pdg = OSH_PART_PDG_NEUTRON;

    /* set the mass of the particle */
    if (!osh_particle_mass_from_pdg(pdg, &p->mass)) {
        p->mass = 0.0;
        p->pdg = OSH_PART_PDG_INVALID;
        p->charge = 0;
        p->is_nucleus = 0;
        p->z = 0;
        p->a = 0;
        return 0;
    }

    if (osh_particle_pdg_is_ion(pdg)) {
        if (!osh_isotope_from_pdg(&iso, pdg)) {
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

    p->pdg = OSH_PART_PDG_INVALID;
    return 0;
}

int osh_particle_pdg_is_ion(int pdg) {
    if (pdg > OSH_PART_PDG_HIBASE)
        return 1;
    else
        return 0;
}

int osh_particle_name_from_pdg(int pdg, char *name_buf, size_t buf_size) {
    struct isotope iso;

    for (size_t i = 0; i < osh_particle_db_len; ++i) {
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
    return 0; /* not found */
}

int osh_particle_symbol_from_pdg(int pdg, char *symbol_buf, size_t buf_size) {
    struct isotope iso;

    for (size_t i = 0; i < osh_particle_db_len; ++i) {
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

    if (!osh_particle_name_from_pdg(p->pdg, name_buf, sizeof(name_buf)))
        name_buf[0] = '\0';

    printf("Particle: %s\n", name_buf);
    printf("----------------------------------------\n");
    printf("  PDG code: %d\n", p->pdg);
    printf("  Mass: %.6f MeV/c^2\n", p->mass);
    printf("  Relative mass: %.6f u\n", p->mass / OSH_AMU);
    printf("  Charge: %d e\n", p->charge);
    if (p->is_nucleus) {
        printf("  Ion: Z=%u, A=%u\n", p->z, p->a);
    } else {
        printf("  Non-ion\n");
    }
    printf("  Run-time properties:\n");
    printf("    Generation: %u\n", p->gen);
    printf("    Nprim: %u\n", p->nprim);
    printf("    Weight: %f\n", p->weight);
    printf("\n");
}