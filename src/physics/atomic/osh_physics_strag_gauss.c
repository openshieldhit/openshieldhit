#include "physics/atomic/osh_physics_strag_gauss.h"

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
 * The literal below equals sqrt(OSH_K_BETHE × OSH_ELECTRON_MASS_MEV) to five
 * digits; it is kept as a literal so the Gaussian result is unchanged by the
 * more precise electron mass now in const.h (used by the Vavilov/Landau path).
 *
 * Reference: PDG "Passage of Particles through Matter" eq. 34.14.
 */
#define BOHR_SQRT_C 0.396128 /* sqrt(K_bethe × m_e c²) = sqrt(0.156917)  [MeV·cm/√(mol/g)·√(g/cm²)] */

double osh_physics_strag_sigma(double z_eff, double z_over_a, double thickness_gcm2) {
    if (z_eff <= 0.0 || z_over_a <= 0.0 || thickness_gcm2 <= 0.0) {
        return 0.0;
    }

    /*
     * σ = sqrt(C_bohr) × z_eff × sqrt((Z/A) × d)
     *
     * Factoring out sqrt(C_bohr) = BOHR_SQRT_C avoids computing the intermediate
     * variance and keeps the expression in one sqrt call.
     */
    return BOHR_SQRT_C * z_eff * sqrt(z_over_a * thickness_gcm2);
}
