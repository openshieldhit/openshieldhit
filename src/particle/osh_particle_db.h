#ifndef OSH_PARTICLE_DB_H
#define OSH_PARTICLE_DB_H

#include <stddef.h>

struct particle_db_entry {
    double mass_mev;
    const char *name;
    const char *symbol;
    int pdg;
    int charge_e;
};

extern struct particle_db_entry const osh_particle_db[];
extern size_t const osh_particle_db_len;

/*
    @brief Get particle database entry by PDG code. Returns NULL if not found.

    @param pdg PDG code of particle to look up
    @return pointer to particle database entry, or NULL if not found

    @author Niels Bassler
*/
struct particle_db_entry const *osh_particle_db_get(int pdg);

#endif /* OSH_PARTICLE_DB_H */