#ifndef OSH_PARTICLE_H
#define OSH_PARTICLE_H

#include <stddef.h>
#include <stdint.h>

/* minimal runtime particle struct for fast cache lookups */
struct particle {
    double mass;   /* MeV */
    double weight; /* statistical weight */

    uint32_t gen;   /* current generation, 0 = primary */
    uint32_t nprim; /* particle root */

    int pdg;            /* PDG number */
    int16_t charge;     /* charge in units of e */
    uint16_t z;         /* atomic number for ions, 0 for non-ions */
    uint16_t a;         /* mass number for ions, 0 for non-ions */
    uint8_t is_nucleus; /* 1 if particle is nucleus, 0 otherwise. All ions except protons and neutrons
                            are considered nuclei. */
};

/* helper functions */

int osh_particle_from_pdg(struct particle *p, int pdg);

int osh_particle_pdg_is_ion(int pdg);

int osh_particle_name_from_pdg(int pdg, char *name_buf, size_t buf_size);
int osh_particle_symbol_from_pdg(int pdg, char *symbol_buf, size_t buf_size);
int osh_particle_mass_from_pdg(int pdg, double *mass_buf);

void osh_print_particle(struct particle const *p);

#endif /* OSH_PARTICLE_H */