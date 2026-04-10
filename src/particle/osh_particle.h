#ifndef OSH_PARTICLE_H
#define OSH_PARTICLE_H

#include <stddef.h>
#include <stdint.h>

/**
 * @brief Species descriptor for a particle type.
 *
 * @details
 * Cold, immutable after setup.  One entry per species in the particle
 * registry; all histories of the same species share a pointer to the same
 * struct particle.  The struct therefore encodes only what is constant for
 * every instance of that species: rest mass, identity codes, and nuclear
 * charge.
 *
 * Per-history runtime state — statistical weight, generation number, primary
 * ancestor index — is NOT stored here.  It belongs in the particle pool (SoA
 * arrays, one entry per live history) or, at scoring time, in struct step.
 * Keeping these concerns separate means:
 *   - struct particle can be held const and shared across threads without
 *     synchronisation.
 *   - The pool's SoA layout can store per-history scalars contiguously,
 *     enabling SIMD access without interleaving cold species metadata.
 *
 * mass is the bare nuclear rest mass (fully stripped ion, no electrons):
 *   M_nuclear = amass [amu] * OSH_AMU  -  Z * m_electron
 * This is the correct relativistic mass for transport kinematics.  Do NOT
 * use A * 940 MeV (free-nucleon approximation) as it introduces ~1% error.
 * Use osh_particle_nuclear_mass_from_pdg() or
 * osh_particle_nuclear_mass_mev_from_za() to populate this field.
 *
 * is_nucleus is 1 for fully stripped nuclei heavier than the proton and
 * neutron (i.e. Z >= 2 or A > 1 excluding neutrons).  The transport engine
 * uses this flag to select the Hubert effective-charge correction in the
 * Bethe-Bloch formula.
 */
struct particle {
    double mass;        /* nuclear rest mass [MeV/c²]; see note above */
    int pdg;            /* PDG Monte Carlo particle numbering scheme code */
    int16_t charge;     /* electric charge in units of e (signed) */
    uint16_t z;         /* atomic number; 0 for non-ion particles */
    uint16_t a;         /* mass number; 0 for non-ion particles */
    uint8_t is_nucleus; /* 1 if heavy nucleus (Z>=2 or A>1, not neutron); 0 otherwise */
};

/* helper functions */

int osh_particle_from_pdg(struct particle *p, int pdg);
int osh_particle_pdg_from_name(char const *name, int *pdg_out);
int osh_particle_from_name(struct particle *p, char const *name);

int osh_particle_pdg_is_ion(int pdg);

int osh_particle_name_from_pdg(int pdg, char *name_buf, size_t buf_size);
int osh_particle_symbol_from_pdg(int pdg, char *symbol_buf, size_t buf_size);

/**
 * @brief Return the default (most abundant) isotope mass number for an element.
 *
 * @details
 * This is the representative isotope used when a dense projectile list is
 * built from atomic number alone, for example when mapping contiguous
 * LOADDEDX columns or constructing the default runtime ion set.
 *
 * @param[in]  z      Atomic number.
 * @param[out] a_out  Receives the default isotope mass number.
 *
 * @returns 1 on success, 0 if Z is outside the isotope database coverage.
 */
int osh_particle_default_isotope_a(unsigned int z, unsigned int *a_out);

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
