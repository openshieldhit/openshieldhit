/**
 * @file fbu_scan.c
 * @brief Scan the Fermi break-up model over excitation energy for one nuclide.
 *
 * @details
 * Reproduces the observables of the standalone Geant4 FermiBreakUp test by
 * Igor Pshenichnov (FermiTest, 2006; used with his kind permission): a parent
 * nucleus (Z, A) at rest is de-excited at fixed E* and the mean fragment
 * multiplicity and mean per-event charge yields are tabulated against the
 * excitation energy per nucleon (100 bins over 0..10 MeV/nucleon, matching
 * the reference histograms in g4fbu_9.1_fixed/).
 *
 * Usage:  fbu_scan <Z> <A>  > <nuc>_osh.dat
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "particle/osh_particle.h"
#include "physics/nuclear/osh_nuclear_fermi_breakup.h"
#include "physics/nuclear/osh_nuclear_handler.h"
#include "random/osh_rng.h"

#define FBU_SCAN_NBINS 100
#define FBU_SCAN_EMAX_PER_A 10.0
#define FBU_SCAN_EVENTS_PER_BIN 2000
#define FBU_SCAN_ZMAX 8

int main(int argc, char **argv) {
    struct osh_nuclear_fermi_breakup model;
    struct osh_rng rng;
    struct osh_nuclear_event ev;
    double mult_sum[FBU_SCAN_NBINS];
    double zyield[FBU_SCAN_NBINS][FBU_SCAN_ZMAX];
    double e_star;
    unsigned int z_par;
    unsigned int a_par;
    int bin;
    int i;
    int zz;
    size_t j;

    if (argc != 3) {
        fprintf(stderr, "usage: fbu_scan <Z> <A>\n");
        return 2;
    }
    z_par = (unsigned int) atoi(argv[1]);
    a_par = (unsigned int) atoi(argv[2]);

    memset(&model, 0, sizeof(model));
    if (osh_nuclear_fermi_breakup_compile(&model) != OSH_OK) {
        fprintf(stderr, "fbu_scan: model compile failed\n");
        return 1;
    }
    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 4242u, 0u);
    memset(mult_sum, 0, sizeof(mult_sum));
    memset(zyield, 0, sizeof(zyield));

    for (bin = 0; bin < FBU_SCAN_NBINS; ++bin) {
        e_star = (bin + 0.5) * (FBU_SCAN_EMAX_PER_A / FBU_SCAN_NBINS) * (double) a_par;
        for (i = 0; i < FBU_SCAN_EVENTS_PER_BIN; ++i) {
            memset(&ev, 0, sizeof(ev));
            ev.kind = OSH_NUCLEAR_EVENT_ABRASION;
            ev.n_fragments = 1u;
            ev.fragments[0].z = z_par;
            ev.fragments[0].a = a_par;
            ev.fragments[0].excitation_energy = e_star;
            osh_nuclear_fermi_breakup_step(&model, &ev.fragments[0], &rng, &ev);

            /* G4 FermiTest counts every final fragment, including an
             * un-broken parent (multiplicity 1 below threshold). */
            mult_sum[bin] += (double) (ev.n_secondaries + ev.n_fragments);
            for (j = 0u; j < ev.n_secondaries; ++j) {
                zz = (ev.secondaries[j].species->pdg == 2212) ? 1 : (int) ev.secondaries[j].species->z;
                if (zz >= 0 && zz < FBU_SCAN_ZMAX) {
                    zyield[bin][zz] += 1.0;
                }
            }
            for (j = 0u; j < ev.n_fragments; ++j) {
                zz = (int) ev.fragments[j].z;
                if (zz >= 0 && zz < FBU_SCAN_ZMAX) {
                    zyield[bin][zz] += 1.0;
                }
            }
        }
    }

    printf("# openshieldhit statistical Fermi break-up — Z=%u A=%u at rest\n", z_par, a_par);
    printf("# E*_per_A_MeV  mean_mult  yield(Z=0..%d)\n", FBU_SCAN_ZMAX - 1);
    for (bin = 0; bin < FBU_SCAN_NBINS; ++bin) {
        printf("%8.4f %10.5f",
               (bin + 0.5) * (FBU_SCAN_EMAX_PER_A / FBU_SCAN_NBINS),
               mult_sum[bin] / FBU_SCAN_EVENTS_PER_BIN);
        for (i = 0; i < FBU_SCAN_ZMAX; ++i) {
            printf(" %10.5f", zyield[bin][i] / FBU_SCAN_EVENTS_PER_BIN);
        }
        printf("\n");
    }

    osh_nuclear_fermi_breakup_free(&model);
    return 0;
}
