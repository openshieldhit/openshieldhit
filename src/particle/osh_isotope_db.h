#ifndef OSH_ISOTOPE_DB_H
#define OSH_ISOTOPE_DB_H

#include "particle/osh_isotope_db_generated.h"

int osh_isotope_from_za(struct isotope *iso, unsigned int z, unsigned a);
int osh_isotope_from_pdg(struct isotope *iso, int pdg);

#endif /* OSH_ISOTOPE_DB_H */