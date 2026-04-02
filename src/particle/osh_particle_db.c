#include "particle/osh_particle_db.h"

#include <stddef.h>

#include "particle/osh_particle_const.h"
#include "particle/osh_particle_pdg.h"

#define MEV(u_) ((double) (u_) * OSH_AMU)

/* Particle Masses, for ions/nucleons these are nuclear, not atomic.
 * - most entries specified in u and converted to MeV at compile time
 * - light nuclei specified directly in MeV as nuclear masses.
 */
struct particle_db_entry const osh_particle_db[] = {
    {OSH_PART_MASS_NEUTRON, "neutron", "n", OSH_PART_PDG_NEUTRON, 0},
    {OSH_PART_MASS_PROTON, "proton", "p", OSH_PART_PDG_PROTON, 1},
    {MEV(1.498347433898e-01), "pi+", "pi+", OSH_PART_PDG_PIPLUS, 1},
    {MEV(1.498347433898e-01), "pi-", "pi-", OSH_PART_PDG_PIMINUS, -1},
    {MEV(1.449033326791e-01), "pi0", "pi0", OSH_PART_PDG_PIZERO, 0},
    {OSH_PART_MASS_NEUTRON, "antineutron", "n~", OSH_PART_PDG_ANEUTRON, 0},
    {OSH_PART_MASS_PROTON, "antiproton", "p~", OSH_PART_PDG_APROTON, -1},

    {MEV(5.299732956210e-01), "K+", "K+", OSH_PART_PDG_KPLUS, 1},
    {MEV(5.299732956210e-01), "K-", "K-", OSH_PART_PDG_KMINUS, -1},

    {MEV(5.342073535546e-01), "K0", "K0", OSH_PART_PDG_KZERO, 0},
    {MEV(5.342073535546e-01), "K~", "K~", OSH_PART_PDG_KTLD, 0},
    {0.0, "gamma", "gamma", OSH_PART_PDG_GAMMA, 0},
    {OSH_PART_MASS_ELECTRON, "electron", "e-", OSH_PART_PDG_ELECTRON, -1},
    {OSH_PART_MASS_ELECTRON, "positron", "e+", OSH_PART_PDG_POSITRON, 1},
    {MEV(1.134289246470e-01), "mu-", "mu-", OSH_PART_PDG_MUMINUS, -1},
    {MEV(1.134289246470e-01), "mu+", "mu+", OSH_PART_PDG_MUPLUS, 1},

    {0.0, "e-neutrino", "nu_e", OSH_PART_PDG_ENU, 0},
    {0.0, "e-antineutrino", "~nu_e", OSH_PART_PDG_EANU, 0},
    {0.0, "mu-neutrino", "nu_mu", OSH_PART_PDG_MUNU, 0},
    {0.0, "mu-antineutrino", "~nu_mu", OSH_PART_PDG_MUANU, 0},

    {1875.61294257, "deuteron", "d", OSH_PART_PDG_DEUTERON, 1},
    {2808.92113298, "triton", "t", OSH_PART_PDG_TRITON, 1},
    {2808.39160743, "He-3", "He3", OSH_PART_PDG_HE3, 2},
    {3727.379378, "He-4", "He4", OSH_PART_PDG_HE4, 2},

    {0.0, "none", "none", OSH_PART_PDG_NONE, 0},
    {0.0, "invalid", "invalid", OSH_PART_PDG_INVALID, 0}

};

size_t const osh_particle_db_len = sizeof(osh_particle_db) / sizeof(osh_particle_db[0]);
