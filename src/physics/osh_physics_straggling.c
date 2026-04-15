#include "physics/osh_physics_straggling.h"

#include <math.h>

/*
 * C_bohr = K_bethe × m_e c²
 *        = 0.307075 [MeV cm²/mol] × 0.511 [MeV]
 *        = 0.156917 [MeV² cm²/mol]
 *
 * Bohr variance:  σ² = C_bohr × z_eff² × (Z/A) × d   [MeV²]
 *                 σ  = sqrt(C_bohr) × z_eff × sqrt((Z/A) × d)
 *
 * where d = ρ·ds is the areal density of the step [g/cm²] and Z/A is in
 * [mol/g], making σ² dimensionally consistent.
 *
 * Reference: PDG "Passage of Particles through Matter" eq. 34.14.
 */
#define BOHR_C 0.156917 /* MeV² cm²/mol */
#define BOHR_SQRT_C 0.396128 /* sqrt(0.156917) */

double osh_physics_straggling_sigma(double z_eff, double z_over_a, double thickness_gcm2) {
    double variance;

    if (z_eff <= 0.0 || z_over_a <= 0.0 || thickness_gcm2 <= 0.0) {
        return 0.0;
    }

    /*
     * σ² = C_bohr × z_eff² × (Z/A) × d
     * σ  = sqrt(C_bohr) × |z_eff| × sqrt((Z/A) × d)
     *
     * Written as a single sqrt to avoid two square-root calls.
     */
    variance = BOHR_C * z_eff * z_eff * z_over_a * thickness_gcm2;
    return sqrt(variance);
}
