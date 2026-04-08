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
int osh_particle_pdg_from_name(char const *name, int *pdg_out);
int osh_particle_from_name(struct particle *p, char const *name);

int osh_particle_pdg_is_ion(int pdg);

int osh_particle_name_from_pdg(int pdg, char *name_buf, size_t buf_size);
int osh_particle_symbol_from_pdg(int pdg, char *symbol_buf, size_t buf_size);

/**
 * @brief Return the nuclear (fully-stripped) rest mass of a particle [MeV/c²].
 *
 * @details
 * For ions the nuclear mass is derived as:
 *   M_nuclear = amass [amu] * OSH_AMU  -  Z * m_electron
 * where amass is the atomic mass from the isotope database.
 * This is the mass of the bare nucleus with all electrons removed, which is
 * the physically correct rest mass for a fully-ionized projectile in transport.
 *
 * For non-ion particles (electrons, muons, etc.) the PDG table mass is returned
 * directly.
 *
 * @param[in]  pdg       PDG code.
 * @param[out] mass_out  Receives nuclear rest mass [MeV/c²].
 *
 * @returns 1 on success, 0 if the PDG code is not found.
 */
int osh_particle_nuclear_mass_from_pdg(int pdg, double *mass_out);

/**
 * @brief Return the atomic mass of a nuclide [amu].
 *
 * @details
 * Returns the tabulated atomic mass (neutral atom, including electron masses
 * and binding energies) in atomic mass units [Da].  This is NOT suitable as
 * the projectile rest mass in kinematics — use osh_particle_nuclear_mass_from_za()
 * for that.  Atomic masses are useful for e.g. material Z/A calculations.
 *
 * @param[in]  z         Atomic number.
 * @param[in]  a         Mass number; 0 means natural element (most abundant isotope).
 * @param[out] amass_out Receives atomic mass [amu].
 *
 * @returns 1 on success, 0 if (Z, A) is not in the isotope database.
 */
int osh_particle_atomic_mass_amu_from_za(unsigned int z, unsigned int a, double *amass_out);

/**
 * @brief Return the nuclear (fully-stripped) rest mass of an ion [MeV/c²].
 *
 * @details
 * Equivalent to osh_particle_nuclear_mass_from_pdg() but takes (Z, A) directly,
 * which is more convenient when no PDG code is available (e.g. when constructing
 * a Bethe projectile descriptor from a material element list).
 *
 * Derived as: M = amass [amu] * OSH_AMU - Z * m_electron.
 *
 * @param[in]  z         Atomic number.
 * @param[in]  a         Mass number.
 * @param[out] mass_out  Receives nuclear rest mass [MeV/c²].
 *
 * @returns 1 on success, 0 if (Z, A) is not in the isotope database.
 */
int osh_particle_nuclear_mass_mev_from_za(unsigned int z, unsigned int a, double *mass_out);

void osh_print_particle(struct particle const *p);

#endif /* OSH_PARTICLE_H */
