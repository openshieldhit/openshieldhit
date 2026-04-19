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

static int populate_light_ion_fields(struct particle *p, int pdg) {
    if (!p) {
        return 0;
    }

    switch (pdg) {
    case OSH_PART_PDG_PROTON:
        p->z = 1u;
        p->a = 1u;
        p->is_nucleus = 0u;
        return 1;
    case OSH_PART_PDG_APROTON:
        p->z = 1u;
        p->a = 1u;
        p->is_nucleus = 0u;
        return 1;
    case OSH_PART_PDG_NEUTRON:
        p->z = 0u;
        p->a = 1u;
        p->is_nucleus = 0u;
        return 1;
    case OSH_PART_PDG_ANEUTRON:
        p->z = 0u;
        p->a = 1u;
        p->is_nucleus = 0u;
        return 1;
    case OSH_PART_PDG_DEUTERON:
        p->z = 1u;
        p->a = 2u;
        p->is_nucleus = 1u;
        return 1;
    case OSH_PART_PDG_TRITON:
        p->z = 1u;
        p->a = 3u;
        p->is_nucleus = 1u;
        return 1;
    case OSH_PART_PDG_HE3:
        p->z = 2u;
        p->a = 3u;
        p->is_nucleus = 1u;
        return 1;
    case OSH_PART_PDG_HE4:
        p->z = 2u;
        p->a = 4u;
        p->is_nucleus = 1u;
        return 1;
    default:
        return 0;
    }
}

int osh_particle_from_pdg(struct particle *p, int pdg) {

    size_t i;
    struct isotope iso;

    if (pdg == OSH_PART_PDG_PROTON_ION) { /* hypothetical hydrogen nucleus */
        pdg = OSH_PART_PDG_PROTON;
    }

    if (pdg == OSH_PART_PDG_NEUTRON_ION) { /* hypothetical neutron nucleus */
        pdg = OSH_PART_PDG_NEUTRON;
    }

    /* set the mass of the particle */
    if (!osh_particle_nuclear_mass_from_pdg(pdg, &p->mass)) {
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
                (void) populate_light_ion_fields(p, pdg);
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

int osh_particle_default_isotope_a(unsigned int z, unsigned int *a_out) {
    unsigned int idx;

    if (!a_out || z >= OSH_ISOTOPE_DB_NELEM) {
        return 0;
    }

    idx = osh_isotopes_idx_default[z];
    if (idx == OSH_ISOTOPE_DB_ERR || idx >= OSH_ISOTOPE_DB_NISO) {
        return 0;
    }

    *a_out = osh_isotope_db[idx].a;
    return 1;
}

int osh_particle_nuclear_mass_from_pdg(int pdg, double *mass_out) {
    size_t i;
    struct isotope iso;

    /* Non-ion particles: use PDG table mass directly. */
    for (i = 0; i < osh_particle_db_len; ++i) {
        if (osh_particle_db[i].pdg == pdg) {
            *mass_out = osh_particle_db[i].mass_mev;
            return 1;
        }
    }

    /*
     * Ions: derive nuclear (fully-stripped) mass from the isotope database.
     *
     *   M_nuclear = amass [amu] * OSH_AMU  -  Z * m_electron
     *
     * The PDG table is checked first (above) for light ions (Z=1,2) where
     * more precise measured masses are available and electron-binding
     * corrections matter at the 10^-4 level.
     */
    if (osh_particle_pdg_is_ion(pdg)) {
        if (osh_isotope_from_pdg(&iso, pdg)) {
            *mass_out = iso.amass * OSH_AMU - (double) iso.z * OSH_PART_MASS_ELECTRON;
            return 1;
        }
    }

    return 0;
}

int osh_particle_atomic_mass_amu_from_za(unsigned int z, unsigned int a, double *amass_out) {
    struct isotope iso;

    if (!amass_out)
        return 0;

    if (a == 0u) {
        /* Natural element: use the default (most abundant) isotope. */
        unsigned int idx = osh_isotopes_idx_default[z];
        if (z >= OSH_ISOTOPE_DB_NELEM || idx == OSH_ISOTOPE_DB_ERR)
            return 0;
        *amass_out = osh_isotope_db[idx].amass;
        return 1;
    }

    if (!osh_isotope_from_za(&iso, z, a))
        return 0;
    *amass_out = iso.amass;
    return 1;
}

int osh_particle_nuclear_mass_mev_from_za(unsigned int z, unsigned int a, double *mass_out) {
    double amass;

    if (!mass_out)
        return 0;
    if (!osh_particle_atomic_mass_amu_from_za(z, a, &amass))
        return 0;

    *mass_out = amass * OSH_AMU - (double) z * OSH_PART_MASS_ELECTRON;
    return 1;
}

void osh_print_particle(struct particle const *p, struct osh_diag_sink const *diag) {
    char name_buf[64];

    if (!p || !diag || !diag->emit || diag->min_level > OSH_DIAG_LEVEL_INFO) {
        return;
    }
    if (!osh_particle_name_from_pdg(p->pdg, name_buf, sizeof(name_buf))) {
        name_buf[0] = '\0';
    }

    OSH_DIAG_INFOF(diag, "Particle: %s", name_buf);
    OSH_DIAG_INFOF(diag, "%s", "------------------------------------------------------------");
    OSH_DIAG_INFOF(diag, "%-18s : %i", "PDG code", p->pdg);
    OSH_DIAG_INFOF(diag, "%-18s : %i", "Z", (int) p->z);
    OSH_DIAG_INFOF(diag, "%-18s : %i", "A", (int) p->a);
    /* Mass is the bare nuclear mass (fully stripped ion, no electrons).
     * Derived as: atomic_mass * OSH_AMU - Z * m_electron (CODATA 2018). */
    OSH_DIAG_INFOF(diag, "%-18s : %-12.5f MeV/c^2  (nuclear, CODATA 2018)", "Mass", p->mass);
    OSH_DIAG_INFOF(diag, "%-18s : %-12.5f Da        (nuclear, CODATA 2018)", "Mass", p->mass / OSH_AMU);
    OSH_DIAG_INFOF(diag, "%-18s : %i e", "Charge", (int) p->charge);
}
