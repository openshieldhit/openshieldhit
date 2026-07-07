/**
 * @file nucre_scan.c
 * @brief Production-spectrum scan of the proton inelastic event generator.
 *
 * @details
 * Stage-0 instrument of the fast nuclear reaction stage work (issues #221 /
 * #260): drives the inelastic branch of the nuclear handler — abrasion
 * followed by Fermi break-up de-excitation, exactly as wired in
 * osh_nuclear_handler_step() — for p + (Z, A) at a fixed incident energy,
 * without any transport.  Tabulated at emission ("production" observables):
 *
 *   - per-species (n, p, d, t, He-3, He-4) yields per inelastic event,
 *     mean kinetic energies, and kinetic-energy histograms;
 *   - the prefragment excitation-energy distribution BEFORE de-excitation
 *     (the E-star supply that feeds the break-up stage);
 *   - leftover unprocessed fragments and their residual excitation after
 *     break-up (energy the transport layer point-discards today).
 *
 * The histogram axis matches the NUCRE reference decks (150 log bins,
 * 0.1..300 MeV) so the curves overlay directly with the plateau spectra of
 * tools/plot_nucre.py.
 *
 * Usage:  nucre_scan <Z> <A> <T_MeV> [nevents] [seed]  > out.dat
 */

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "particle/osh_particle.h"
#include "particle/osh_particle_pdg.h"
#include "physics/nuclear/osh_nuclear_abrasion.h"
#include "physics/nuclear/osh_nuclear_fermi_breakup.h"
#include "physics/nuclear/osh_nuclear_handler.h"
#include "physics/nuclear/osh_nuclear_tripathi.h"
#include "random/osh_rng.h"

#define NUCRE_SCAN_NBINS 150
#define NUCRE_SCAN_ELO_MEV 0.1
#define NUCRE_SCAN_EHI_MEV 300.0
#define NUCRE_SCAN_NSPECIES 6
#define NUCRE_SCAN_DEFAULT_EVENTS 200000UL
#define NUCRE_SCAN_DEFAULT_SEED 4242u

static char const *const s_species_name[NUCRE_SCAN_NSPECIES] = {
    "n",
    "p",
    "d",
    "t",
    "he3",
    "alpha",
};

/** @brief Parse a base-10 unsigned integer argument without sign wrapping. */
static int parse_uint_arg(char const *text, unsigned int *out) {
    char *end;
    unsigned long value;

    if (!text || !out || text[0] == '-' || text[0] == '+') {
        return 0;
    }
    errno = 0;
    end = NULL;
    value = strtoul(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || value > (unsigned long) UINT_MAX) {
        return 0;
    }
    *out = (unsigned int) value;
    return 1;
}

/** @brief Parse a base-10 unsigned long argument without sign wrapping. */
static int parse_ulong_arg(char const *text, unsigned long *out) {
    char *end;
    unsigned long value;

    if (!text || !out || text[0] == '-' || text[0] == '+') {
        return 0;
    }
    errno = 0;
    end = NULL;
    value = strtoul(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0') {
        return 0;
    }
    *out = value;
    return 1;
}

/** @brief Parse a positive finite floating-point argument. */
static int parse_positive_double_arg(char const *text, double *out) {
    char *end;
    double value;

    if (!text || !out) {
        return 0;
    }
    errno = 0;
    end = NULL;
    value = strtod(text, &end);
    if (errno != 0 || end == text || *end != '\0' || !isfinite(value) || !(value > 0.0)) {
        return 0;
    }
    *out = value;
    return 1;
}

/* Map a secondary's species to the fixed tally row; -1 for non-whitelist. */
static int species_index(struct particle const *sp) {
    if (sp->pdg == OSH_PART_PDG_NEUTRON) {
        return 0;
    }
    if (sp->pdg == OSH_PART_PDG_PROTON) {
        return 1;
    }
    if (sp->z == 1u && sp->a == 2u) {
        return 2;
    }
    if (sp->z == 1u && sp->a == 3u) {
        return 3;
    }
    if (sp->z == 2u && sp->a == 3u) {
        return 4;
    }
    if (sp->z == 2u && sp->a == 4u) {
        return 5;
    }
    return -1;
}

/* Log-axis bin index, or -1 below range / NUCRE_SCAN_NBINS overflow clamp. */
static int energy_bin(double e_mev) {
    double u;
    if (e_mev <= NUCRE_SCAN_ELO_MEV) {
        return -1;
    }
    u = log(e_mev / NUCRE_SCAN_ELO_MEV) / log(NUCRE_SCAN_EHI_MEV / NUCRE_SCAN_ELO_MEV);
    if (u >= 1.0) {
        return NUCRE_SCAN_NBINS - 1;
    }
    return (int) (u * NUCRE_SCAN_NBINS);
}

int main(int argc, char **argv) {
    struct osh_nuclear_fermi_breakup fbu;
    struct osh_rng rng;
    struct osh_nuclear_event ev;
    static double hist[NUCRE_SCAN_NSPECIES][NUCRE_SCAN_NBINS];
    static double estar_hist[NUCRE_SCAN_NBINS];
    double yield[NUCRE_SCAN_NSPECIES];
    double esum[NUCRE_SCAN_NSPECIES];
    double below[NUCRE_SCAN_NSPECIES];
    double dir[3];
    double sigma_cm2;
    double t_mev;
    double estar;
    double estar_sum;
    double estar_sq_sum;
    double estar_max;
    double leftover_frag;
    double leftover_estar;
    double edges_ratio;
    double elo;
    unsigned long nevents;
    unsigned long n_with_fragment;
    unsigned long iev;
    unsigned int z_tgt;
    unsigned int a_tgt;
    unsigned int seed;
    size_t j;
    int s;
    int bin;

    if (argc < 4 || argc > 6) {
        fprintf(stderr, "usage: nucre_scan <Z> <A> <T_MeV> [nevents] [seed]\n");
        return 2;
    }
    if (!parse_uint_arg(argv[1], &z_tgt) || !parse_uint_arg(argv[2], &a_tgt)
        || !parse_positive_double_arg(argv[3], &t_mev)) {
        fprintf(stderr, "nucre_scan: invalid target/energy/event count\n");
        return 2;
    }
    nevents = NUCRE_SCAN_DEFAULT_EVENTS;
    seed = NUCRE_SCAN_DEFAULT_SEED;
    if (argc > 4 && !parse_ulong_arg(argv[4], &nevents)) {
        fprintf(stderr, "nucre_scan: invalid target/energy/event count\n");
        return 2;
    }
    if (argc > 5 && !parse_uint_arg(argv[5], &seed)) {
        fprintf(stderr, "nucre_scan: invalid target/energy/event count\n");
        return 2;
    }
    if (z_tgt == 0u || a_tgt < 2u || nevents == 0UL) {
        fprintf(stderr, "nucre_scan: invalid target/energy/event count\n");
        return 2;
    }

    sigma_cm2 = osh_nuclear_tripathi_sigma(1u, 1u, (double) z_tgt, (double) a_tgt, t_mev);
    if (sigma_cm2 <= 0.0) {
        fprintf(stderr, "nucre_scan: Tripathi sigma is zero (below Coulomb threshold?)\n");
        return 1;
    }

    memset(&fbu, 0, sizeof(fbu));
    if (osh_nuclear_fermi_breakup_compile(&fbu) != OSH_OK) {
        fprintf(stderr, "nucre_scan: Fermi break-up compile failed\n");
        return 1;
    }
    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, seed, 0u);

    memset(hist, 0, sizeof(hist));
    memset(estar_hist, 0, sizeof(estar_hist));
    memset(yield, 0, sizeof(yield));
    memset(esum, 0, sizeof(esum));
    memset(below, 0, sizeof(below));
    dir[0] = 0.0;
    dir[1] = 0.0;
    dir[2] = 1.0;
    estar_sum = 0.0;
    estar_sq_sum = 0.0;
    estar_max = 0.0;
    leftover_frag = 0.0;
    leftover_estar = 0.0;
    n_with_fragment = 0UL;

    for (iev = 0UL; iev < nevents; ++iev) {
        memset(&ev, 0, sizeof(ev));
        osh_nuclear_abrasion_step(t_mev, dir, (double) a_tgt, (double) z_tgt, sigma_cm2, &rng, &ev);

        /* Prefragment excitation BEFORE de-excitation: the E-star supply. */
        if (ev.n_fragments > 0u) {
            estar = ev.fragments[0].excitation_energy;
            ++n_with_fragment;
            estar_sum += estar;
            estar_sq_sum += estar * estar;
            if (estar > estar_max) {
                estar_max = estar;
            }
            bin = energy_bin(estar);
            if (bin >= 0) {
                estar_hist[bin] += 1.0;
            }
            osh_nuclear_fermi_breakup_step(&fbu, &ev.fragments[0], &rng, &ev);
        }

        /* Secondaries at emission (abrasion knockouts, the escaping cascade
         * proton, and break-up products). */
        for (j = 0u; j < ev.n_secondaries; ++j) {
            s = species_index(ev.secondaries[j].species);
            if (s < 0) {
                continue;
            }
            yield[s] += 1.0;
            esum[s] += ev.secondaries[j].energy;
            bin = energy_bin(ev.secondaries[j].energy);
            if (bin >= 0) {
                hist[s][bin] += 1.0;
            } else {
                below[s] += 1.0;
            }
        }

        /* Unprocessed residues after break-up; their excitation is what the
         * transport layer currently discards on pool injection. */
        leftover_frag += (double) ev.n_fragments;
        for (j = 0u; j < ev.n_fragments; ++j) {
            leftover_estar += ev.fragments[j].excitation_energy;
        }
    }

    printf("# nucre_scan — proton inelastic production spectra (abrasion + Fermi break-up)\n");
    printf("# target Z=%u A=%u  T_lab=%.6g MeV  sigma_inel=%.3f mb  events=%lu  seed=%u\n",
           z_tgt,
           a_tgt,
           t_mev,
           sigma_cm2 * 1.0e27,
           nevents,
           seed);
    printf("# species  yield/event  mean_E_MeV  frac_below_%.3gMeV\n", NUCRE_SCAN_ELO_MEV);
    for (s = 0; s < NUCRE_SCAN_NSPECIES; ++s) {
        printf("#   %-5s  %10.5f  %10.4f  %10.5f\n",
               s_species_name[s],
               yield[s] / (double) nevents,
               (yield[s] > 0.0) ? esum[s] / yield[s] : 0.0,
               (yield[s] > 0.0) ? below[s] / yield[s] : 0.0);
    }
    printf("# prefragment: frac_with_fragment=%.5f  Estar_mean=%.4f  Estar_std=%.4f  Estar_max=%.2f MeV\n",
           (double) n_with_fragment / (double) nevents,
           (n_with_fragment > 0UL) ? estar_sum / (double) n_with_fragment : 0.0,
           (n_with_fragment > 0UL)
               ? sqrt(fmax(0.0,
                           estar_sq_sum / (double) n_with_fragment
                               - (estar_sum / (double) n_with_fragment) * (estar_sum / (double) n_with_fragment)))
               : 0.0,
           estar_max);
    printf("# after break-up: leftover_fragments/event=%.5f  leftover_Estar/event=%.5f MeV (discarded by transport)\n",
           leftover_frag / (double) nevents,
           leftover_estar / (double) nevents);
    printf("# columns: Elo_MeV Ehi_MeV n p d t he3 alpha prefrag_Estar   (counts/event per bin)\n");

    edges_ratio = pow(NUCRE_SCAN_EHI_MEV / NUCRE_SCAN_ELO_MEV, 1.0 / NUCRE_SCAN_NBINS);
    for (bin = 0; bin < NUCRE_SCAN_NBINS; ++bin) {
        elo = NUCRE_SCAN_ELO_MEV * pow(edges_ratio, (double) bin);
        printf("%12.6g %12.6g", elo, elo * edges_ratio);
        for (s = 0; s < NUCRE_SCAN_NSPECIES; ++s) {
            printf(" %12.6g", hist[s][bin] / (double) nevents);
        }
        printf(" %12.6g\n", estar_hist[bin] / (double) nevents);
    }

    osh_nuclear_fermi_breakup_free(&fbu);
    return 0;
}
