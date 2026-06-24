#include "physics/nuclear/osh_nuclear_compound.h"

#include <string.h>

#include "common/osh_diag.h"
#include "physics/nuclear/osh_nuclear_fermi_breakup.h"
#include "physics/nuclear/osh_nuclear_handler.h"

/*
 * One-time heavy-A sink warning table.
 * Not thread-safe; acceptable for the current single-threaded transport model.
 * This is the future SMM entry point — the sink is intentionally minimal.
 */
#define MAX_HEAVY_WARNED 32
static unsigned int s_warned_z[MAX_HEAVY_WARNED];
static unsigned int s_warned_a[MAX_HEAVY_WARNED];
static int s_n_warned = 0;

static void warn_heavy_once(unsigned int z, unsigned int a,
                            struct osh_diag_sink const *diag) {
    int i;
    int found;

    if (s_n_warned >= MAX_HEAVY_WARNED)
        return;

    found = 0;
    for (i = 0; i < s_n_warned; ++i) {
        if (s_warned_z[i] == z && s_warned_a[i] == a) {
            found = 1;
            break;
        }
    }
    if (found)
        return;

    OSH_DIAG_INFOF(diag,
        "compound nucleus: A=%u Z=%u exceeds FBU domain (A > %d); "
        "excitation energy deposited locally (heavy-A sink, future SMM)",
        a, z, OSH_FERMI_BREAKUP_AMAX);

    s_warned_z[s_n_warned] = z;
    s_warned_a[s_n_warned] = a;
    s_n_warned++;
}

void osh_nuclear_compound_step(unsigned int z, unsigned int a,
                               double e_star_mev, double const p_lab_mev[3],
                               struct osh_nuclear_fermi_breakup const *fbu,
                               struct osh_diag_sink const *diag,
                               struct osh_rng *rng,
                               struct osh_nuclear_event *event_out) {
    struct osh_nuclear_fragment fragment;

    memset(event_out, 0, sizeof(*event_out));

    /* -- light compound nucleus: route to Fermi break-up ------------------- */
    if (a <= (unsigned int)OSH_FERMI_BREAKUP_AMAX) {
        fragment.z                = z;
        fragment.a                = a;
        fragment.excitation_energy = e_star_mev;
        fragment.p[0]             = p_lab_mev[0];
        fragment.p[1]             = p_lab_mev[1];
        fragment.p[2]             = p_lab_mev[2];
        osh_nuclear_fermi_breakup_step(fbu, &fragment, rng, event_out);
        return;
    }

    /* -- heavy-A sink: deposit excitation locally (future SMM seam) -------- */
    warn_heavy_once(z, a, diag);
    event_out->kind        = OSH_NUCLEAR_EVENT_ABSORB;
    event_out->n_secondaries = 0u;
    event_out->n_fragments   = 0u;
}
