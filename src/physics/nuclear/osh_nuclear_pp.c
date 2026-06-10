#include "physics/nuclear/osh_nuclear_pp.h"

#include "common/osh_interpolate.h"
#include "physics/nuclear/osh_nuclear_tripathi.h"
#include "random/osh_rng.h"

#include "physics/nuclear/osh_nuclear_pp_data.h"

/* ---- Cross section ------------------------------------------------------- */

double osh_nuclear_pp_sigma_el(double e_lab_mev) {
    float e_f;
    long int ei;
    float te;
    float sigma_mb;

    if (e_lab_mev < (double) osh_pp_elab_mev[0]) {
        return 0.0;
    }

    e_f = (float) e_lab_mev;
    ei  = osh_binary_search_f(e_f, osh_pp_elab_mev, OSH_PP_NENERGY);
    if (ei < 0) {
        return 0.0;
    }

    te = (osh_pp_elab_mev[ei + 1] > osh_pp_elab_mev[ei])
           ? (e_f - osh_pp_elab_mev[ei]) / (osh_pp_elab_mev[ei + 1] - osh_pp_elab_mev[ei])
           : 0.0f;
    if (te < 0.0f) te = 0.0f;
    if (te > 1.0f) te = 1.0f;

    sigma_mb = (1.0f - te) * osh_pp_sigma_tot_mb[ei] * (1.0f - osh_pp_sigma_in_frac[ei])
             + te * osh_pp_sigma_tot_mb[ei + 1] * (1.0f - osh_pp_sigma_in_frac[ei + 1]);

    return (double) sigma_mb * 1.0e-27; /* mb → cm² */
}

/* ---- CDF angle sampling -------------------------------------------------- */

double osh_nuclear_pp_sample_cos_theta_cm(double e_lab_mev, struct osh_rng *rng) {
    float e_f;
    long int ei;
    float te;
    double U;
    int j;
    float cdf_next;
    float cdf_lo, cdf_hi;
    float width, ta;

    e_f = (float) e_lab_mev;
    ei  = osh_binary_search_f(e_f, osh_pp_elab_mev, OSH_PP_NENERGY);
    if (ei < 0) {
        ei = 0;
    }

    te = (osh_pp_elab_mev[ei + 1] > osh_pp_elab_mev[ei])
           ? (e_f - osh_pp_elab_mev[ei]) / (osh_pp_elab_mev[ei + 1] - osh_pp_elab_mev[ei])
           : 0.0f;
    if (te < 0.0f) te = 0.0f;
    if (te > 1.0f) te = 1.0f;

    U = osh_rng_double(rng); /* one deviate in [0, 1) */

    /*
     * Forward scan for angle bracket.  osh_binary_search_f is not used here
     * because CDF rows may contain flat segments (equal consecutive entries)
     * which violate its strict-monotonicity precondition.
     * Loop finds lowest j such that interpolated CDF at j+1 >= U.
     * Loop bound is NANGLE-2 so j+1 is always a valid index.
     */
    j = 0;
    while (j < OSH_PP_NANGLE - 2) {
        cdf_next = (1.0f - te) * osh_pp_cdf[ei][j + 1] + te * osh_pp_cdf[ei + 1][j + 1];
        if ((float) U <= cdf_next) {
            break;
        }
        ++j;
    }

    cdf_lo = (1.0f - te) * osh_pp_cdf[ei][j]     + te * osh_pp_cdf[ei + 1][j];
    cdf_hi = (1.0f - te) * osh_pp_cdf[ei][j + 1] + te * osh_pp_cdf[ei + 1][j + 1];

    /* Sub-bin linear interpolation in cos(θ) */
    width = cdf_hi - cdf_lo;
    ta    = (width > 0.0f) ? (float) ((U - (double) cdf_lo) / (double) width) : 0.0f;

    return (double) ((1.0f - ta) * osh_pp_cos_theta_axis[j] + ta * osh_pp_cos_theta_axis[j + 1]);
}

/* ---- Mean free path ------------------------------------------------------ */

double osh_nuclear_pp_lambda_gcm2(double hydrogen_mass_fraction, double sigma_el_cm2) {
    if (hydrogen_mass_fraction <= 0.0 || sigma_el_cm2 <= 0.0) {
        return 1.0e30;
    }
    /*
     * λ = A_H / (f_H × N_A × σ_el)
     * osh_nuclear_lambda_gcm2(A, σ) = A / (N_A × σ), so divide by f_H.
     */
    return osh_nuclear_lambda_gcm2(1.00794, sigma_el_cm2) / hydrogen_mass_fraction;
}
