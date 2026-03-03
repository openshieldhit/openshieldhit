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

#endif /* OSH_PARTICLE_DB_H */