#include "material/osh_material_atomic_data.h"

#include <stddef.h>

#include "particle/osh_particle.h"

/*
 * Natural-element atomic weights [Da].
 *
 * Source: IUPAC standard atomic weights. These values represent natural elements
 * and must not be confused with isotope masses from particle/.
 */
static double const natural_atomic_masses_da[] = {
    1.00794,     4.002602,  6.941,      9.012182, 10.811,      12.0107, 14.0067,   15.9994, 18.9984032, 20.1797,
    22.98976928, 24.3050,   26.9815386, 28.0855,  30.973762,   32.065,  35.453,    39.948,  39.0983,    40.078,
    44.955912,   47.867,    50.9415,    51.9961,  54.938045,   55.845,  58.933195, 58.6934, 63.546,     65.38,
    69.723,      72.64,     74.92160,   78.96,    79.904,      83.798,  85.4678,   87.62,   88.90585,   91.224,
    92.90638,    95.96,     98.0,       101.07,   102.90550,   106.42,  107.8682,  112.411, 114.818,    118.710,
    121.760,     127.60,    126.90447,  131.293,  132.9054519, 137.327, 138.90547, 140.116, 140.90765,  144.242,
    145.0,       150.36,    151.964,    157.25,   158.92535,   162.500, 164.93032, 167.259, 168.93421,  173.054,
    174.9668,    178.49,    180.94788,  183.84,   186.207,     190.23,  192.217,   195.084, 196.966569, 200.59,
    204.3833,    207.2,     208.98040,  209.0,    210.0,       222.0,   223.0,     226.0,   227.0,      232.03806,
    231.03588,   238.02891, 237.0,      244.0,    243.0,       247.0,   247.0,     251.0,   252.0,      257.0,
    258.0,       259.0,     262.0,      267.0,    268.0,       271.0,   272.0,     270.0,   276.0,      281.0,
    280.0,       285.0};

enum osh_status osh_material_natural_atomic_mass_da(unsigned int z, double *mass_out) {
    size_t idx;

    if (!mass_out || z == 0u
        || z > (unsigned int) (sizeof(natural_atomic_masses_da) / sizeof(natural_atomic_masses_da[0]))) {
        return OSH_EINVAL;
    }

    idx = (size_t) z - 1u;
    *mass_out = natural_atomic_masses_da[idx];
    return OSH_OK;
}

enum osh_status osh_material_atomic_mass_da(unsigned int z, unsigned int a, double *mass_out) {
    if (!mass_out || z == 0u) {
        return OSH_EINVAL;
    }

    if (a == 0u) {
        return osh_material_natural_atomic_mass_da(z, mass_out);
    }

    if (!osh_particle_atomic_mass_amu_from_za(z, a, mass_out)) {
        return OSH_EINVAL;
    }

    return OSH_OK;
}
