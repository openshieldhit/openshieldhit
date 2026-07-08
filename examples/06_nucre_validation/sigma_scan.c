/**
 * @file sigma_scan.c
 * @brief Dump the p+A cross sections the transport actually uses (issue #277).
 *
 * @details
 * Prints, for p + (Z, A) over an incident-energy grid, exactly the quantities
 * the nuclear handler feeds into the removal rate:
 *
 *   - sigma_R    : Tripathi reaction cross section [mb]
 *                  (osh_nuclear_tripathi_sigma());
 *   - sigma_el   : p+A nuclear elastic cross section [mb] as configured
 *                  (osh_nuclear_elastic_sigma(), currently
 *                  OSH_NUCLEAR_ELASTIC_SIGMA_FACTOR * sigma_R);
 *   - B          : diffraction slope of the exp(-B|t|) angular model
 *                  [(GeV/c)^-2] (osh_nuclear_elastic_slope());
 *   - theta_med  : median CM scattering angle of that model [deg], from
 *                  |t|_med = ln(2)/B at the CM momentum for this energy.
 *
 * The output overlays directly with the committed reference tables in
 * tests/reference/xsec/ via tools/plot_xsec_pa.py.
 *
 * Usage:  sigma_scan <Z> <A> [e_min_MeV] [e_max_MeV] [n_steps]  > out.dat
 */

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "particle/osh_particle_const.h"
#include "physics/nuclear/osh_nuclear_elastic.h"
#include "physics/nuclear/osh_nuclear_pp.h"
#include "physics/nuclear/osh_nuclear_tripathi.h"

#define SIGMA_SCAN_E_MIN_MEV 5.0
#define SIGMA_SCAN_E_MAX_MEV 250.0
#define SIGMA_SCAN_N_STEPS 246
#define SIGMA_SCAN_CM2_TO_MB 1.0e27

/** CM momentum [MeV/c] of a proton with lab kinetic energy @p t_mev on a
 * target of mass @p m_target_mev (relativistic two-body). */
static double _p_cm_mev(double t_mev, double m_target_mev) {
    double m_p;
    double p_lab;
    double w_inv;

    m_p = OSH_PART_MASS_PROTON;
    p_lab = sqrt(t_mev * (t_mev + (2.0 * m_p)));
    w_inv = sqrt((m_p * m_p) + (m_target_mev * m_target_mev) + (2.0 * m_target_mev * (t_mev + m_p)));
    return (w_inv > 0.0) ? (p_lab * m_target_mev / w_inv) : 0.0;
}

int main(int argc, char **argv) {
    unsigned int zt;
    unsigned int at;
    double e_min;
    double e_max;
    long n_steps;
    long i;
    double m_target;
    double b_slope;

    if (argc < 3) {
        fprintf(stderr, "usage: %s <Z> <A> [e_min_MeV] [e_max_MeV] [n_steps]\n", argv[0]);
        return EXIT_FAILURE;
    }
    errno = 0;
    zt = (unsigned int) strtoul(argv[1], NULL, 10);
    at = (unsigned int) strtoul(argv[2], NULL, 10);
    e_min = (argc > 3) ? strtod(argv[3], NULL) : SIGMA_SCAN_E_MIN_MEV;
    e_max = (argc > 4) ? strtod(argv[4], NULL) : SIGMA_SCAN_E_MAX_MEV;
    n_steps = (argc > 5) ? strtol(argv[5], NULL, 10) : SIGMA_SCAN_N_STEPS;
    if (errno != 0 || zt == 0u || at == 0u || !(e_min > 0.0) || !(e_max > e_min) || n_steps < 2) {
        fprintf(stderr, "sigma_scan: invalid arguments\n");
        return EXIT_FAILURE;
    }

    /* Crude nuclear mass for the CM conversion; the sub-percent binding-energy
     * effect is irrelevant for a median-angle diagnostic. */
    m_target = (double) at * OSH_AMU;
    b_slope = osh_nuclear_elastic_slope((double) at);

    printf("# sigma_scan: p + (Z=%u, A=%u), OSH transport cross sections (issue #277)\n", zt, at);
    printf("# B = %.6e (MeV/c)^-2 (exp(-B|t|) diffraction slope, constant in E)\n", b_slope);
    printf("# columns: E_MeV  sigma_R_mb  sigma_el_mb  theta_cm_med_deg\n");
    for (i = 0; i < n_steps; ++i) {
        double e;
        double sigma_r_cm2;
        double sigma_el_cm2;
        double p_cm;
        double t_med;
        double cos_med;
        double theta_med;

        e = e_min + ((e_max - e_min) * (double) i / (double) (n_steps - 1));
        if (zt == 1u && at == 1u) {
            /* Hydrogen target: the handler routes this to the pp channel
             * (no p+p inelastic below the pion threshold, no p+A elastic). */
            sigma_r_cm2 = 0.0;
            sigma_el_cm2 = osh_nuclear_pp_sigma_el(e);
        } else {
            sigma_r_cm2 = osh_nuclear_tripathi_sigma(1u, 1u, (double) zt, (double) at, e);
            sigma_el_cm2 = osh_nuclear_elastic_sigma(1u, 1u, (double) zt, (double) at, e);
        }

        p_cm = _p_cm_mev(e, m_target);
        theta_med = 0.0;
        /* The exp(-B|t|) angle model applies to the p+A channel only; for the
         * hydrogen row the pp sampler has its own angular distribution. */
        if (!(zt == 1u && at == 1u) && p_cm > 0.0 && b_slope > 0.0) {
            t_med = log(2.0) / b_slope;
            cos_med = 1.0 - (t_med / (2.0 * p_cm * p_cm));
            if (cos_med < -1.0) {
                cos_med = -1.0;
            }
            theta_med = acos(cos_med) * 180.0 / M_PI;
        }

        printf("%10.4f %12.4f %12.4f %10.4f\n",
               e,
               sigma_r_cm2 * SIGMA_SCAN_CM2_TO_MB,
               sigma_el_cm2 * SIGMA_SCAN_CM2_TO_MB,
               theta_med);
    }
    return EXIT_SUCCESS;
}
