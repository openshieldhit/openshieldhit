#include "physics/nuclear/osh_nuclear_sigma_reac.h"

#include <math.h>
#include <stddef.h>

#include "physics/nuclear/osh_nuclear_sigma_reac_data.h"
#include "physics/nuclear/osh_nuclear_tripathi.h"

#define OSH_SIGMA_REAC_MB_TO_CM2 1.0e-27

#define OSH_SIGMA_REAC_NTABLES                                                                                         \
    ((unsigned int) (sizeof(osh_nuclear_sigma_reac_tables) / sizeof(osh_nuclear_sigma_reac_tables[0])))

/** Table for the (z, a) target, or NULL when the target is not tabulated. */
static struct osh_nuclear_sigma_reac_table const *_sigma_reac_table(unsigned int z, unsigned int a) {
    unsigned int i;

    for (i = 0u; i < OSH_SIGMA_REAC_NTABLES; ++i) {
        if (osh_nuclear_sigma_reac_tables[i].z == z && osh_nuclear_sigma_reac_tables[i].a == a) {
            return &osh_nuclear_sigma_reac_tables[i];
        }
    }
    return NULL;
}

/** Lin-lin interpolation on the uniform table grid; flat clamp above the last
 * point (supported by the Renberg 230 MeV data), 0 below the first point. */
static double _sigma_reac_lookup_cm2(struct osh_nuclear_sigma_reac_table const *table, double e_mev) {
    double x;
    double frac;
    unsigned int idx;
    double sigma_mb;

    x = (e_mev - OSH_NUCLEAR_SIGMA_REAC_E_MIN_MEV) / OSH_NUCLEAR_SIGMA_REAC_E_STEP_MEV;
    if (x <= 0.0) {
        sigma_mb = (double) table->sigma_mb[0];
    } else if (x >= (double) (OSH_NUCLEAR_SIGMA_REAC_NPOINTS - 1)) {
        sigma_mb = (double) table->sigma_mb[OSH_NUCLEAR_SIGMA_REAC_NPOINTS - 1];
    } else {
        idx = (unsigned int) x;
        frac = x - (double) idx;
        sigma_mb = (double) table->sigma_mb[idx] + (frac * ((double) table->sigma_mb[idx + 1u] - (double) table->sigma_mb[idx]));
    }
    return sigma_mb * OSH_SIGMA_REAC_MB_TO_CM2;
}

double osh_nuclear_sigma_reac(unsigned int zp, unsigned int ap, double zt, double at, double e_lab_per_nucleon) {
    struct osh_nuclear_sigma_reac_table const *table;
    double zt_round;
    double at_round;

    if (zp == 1u && ap == 1u && zt > 0.0 && at > 0.0) {
        /* Tables are keyed on integer (Z, A); match only exact nuclides so a
         * compound-average pseudo-nucleus never silently hits a table. */
        zt_round = floor(zt + 0.5);
        at_round = floor(at + 0.5);
        if (zt == zt_round && at == at_round) {
            table = _sigma_reac_table((unsigned int) zt_round, (unsigned int) at_round);
            if (table != NULL) {
                return _sigma_reac_lookup_cm2(table, e_lab_per_nucleon);
            }
        }
    }
    return osh_nuclear_tripathi_sigma(zp, ap, zt, at, e_lab_per_nucleon);
}
