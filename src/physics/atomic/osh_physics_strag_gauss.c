#include "physics/atomic/osh_physics_strag_gauss.h"

#include <math.h>

/*
 * C_bohr = K_bethe × m_e c²
 *        = 0.307075 [MeV cm²/mol] × 0.511 [MeV]
 *        = 0.156917 [MeV² cm²/mol]
 *
 * Bohr variance:  σ² = C_bohr × z_eff² × (Z/A) × d × γ² × (1 − β²/2)   [MeV²]
 *                 σ  = sqrt(C_bohr) × z_eff × sqrt((Z/A) × d) × γ × sqrt(1 − β²/2)
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

double osh_physics_strag_sigma(double z_eff, double z_over_a, double thickness_gcm2, double beta2) {
    if (z_eff <= 0.0 || z_over_a <= 0.0 || thickness_gcm2 <= 0.0 || beta2 <= 0.0 || beta2 >= 1.0) {
        return 0.0;
    }

    /*
     * σ = sqrt(C_bohr) × z_eff × sqrt((Z/A) × d × (1 − β²/2) / (1 − β²))
     *
     * This is algebraically the same as multiplying the classical σ by
     * γ × sqrt(1 − β²/2), but folds the relativistic variance factor into the
     * existing square root.
     */
    return BOHR_SQRT_C * z_eff * sqrt(z_over_a * thickness_gcm2 * (1.0 - 0.5 * beta2) / (1.0 - beta2));
}
