#include "scoring/save/osh_scoring_save_mcpl.h"

#include <math.h>
#include <string.h>

#include "mcpl.h"
#include "scoring/runtime/osh_scoring_mcpl_record.h"
#include "scoring/runtime/osh_scoring_output_runtime.h"
#include "scoring/runtime/osh_scoring_runtime.h"

/* mcpl_add_particle() aborts the process (mcpl.c's default error handler) on a
 * direction vector whose squared norm deviates from 1 by more than 1e-5.
 * struct step documents st->w as always a unit vector, so this should be a
 * no-op in practice; it is cheap insurance against a hard process abort from
 * floating-point drift over a long transport chain, paid once per record at
 * save time (never on the hot path that books mcpl_records). A degenerate
 * (near-zero) vector maps to +Z rather than dividing by ~0. */
static void normalize_direction(double dir[3]) {
    double const norm2 = dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2];
    double inv_norm;

    if (!(norm2 > 1.0e-20)) {
        dir[0] = 0.0;
        dir[1] = 0.0;
        dir[2] = 1.0;
        return;
    }
    inv_norm = 1.0 / sqrt(norm2);
    dir[0] *= inv_norm;
    dir[1] *= inv_norm;
    dir[2] *= inv_norm;
}

enum osh_status osh_scoring_save_mcpl_output(struct osh_scoring_workspace const *ws,
                                             struct osh_scoring_runtime const *rt,
                                             unsigned long long nstat,
                                             size_t output_idx) {
    struct osh_scoring_output_runtime const *out;
    struct osh_scoring_page_runtime const *page;
    struct osh_scoring_accumulator const *acc;
    mcpl_outfile_t of;
    mcpl_particle_t particle;
    size_t count;
    size_t i;

    (void) ws;

    if (!rt || output_idx >= rt->noutputs) {
        return OSH_EINVAL;
    }
    if (nstat == 0ull) {
        return OSH_EINVAL;
    }

    out = &rt->outputs[output_idx];
    /* osh_scoring_compile() already enforces this shape for FileFormat MCPL
     * (Phase 2 validation); re-check rather than trust the caller blindly. */
    if (out->npages != 1u) {
        return OSH_ENOTSUP;
    }
    page = &rt->pages[out->page_indices[0]];
    acc = &page->acc;
    if (page->score_kind != OSH_SCORING_SCORE_MCPL || !acc->mcpl_records || !acc->mcpl_count) {
        return OSH_ESTATE;
    }
    count = *acc->mcpl_count;
    if (count > acc->mcpl_capacity) {
        return OSH_ESTATE; /* would indicate a hot-path bound-check bug, not user error */
    }

    /* mcpl_create_outfile()/mcpl_add_particle() have no error-return path: the
     * vendored library aborts the process (mcpl.c's default error handler) on
     * a failure such as an unwritable path. Accepted upstream MCPL behaviour,
     * unlike the enum osh_status convention the rest of this codebase uses —
     * see src/thirdparty/mcpl/README.md. */
    of = mcpl_create_outfile(out->filename);
    mcpl_hdr_set_srcname(of, "openshieldhit");
    mcpl_enable_doubleprec(of);
    mcpl_enable_userflags(of);
    mcpl_hdr_add_comment(of, "userflags: generation number (0 = beam primary, N = Nth-generation secondary)");
    /* Known in full at save time (unlike a live/streaming writer that has to
     * reserve this and patch it in once the run finishes): the total primary
     * count this run represents, letting a downstream consumer scale a partial
     * phase-space dump back to a per-primary basis.
     *
     * MCPL standardises no key name for this, so every producer picks its own
     * (Geant4 phase-space writers commonly use "launched_primaries"); a
     * consumer reading the wrong key silently mis-scales every downstream
     * result rather than failing.  Name the key in a header comment as well as
     * the stat, so the file says what it means without reference to these docs. */
    mcpl_hdr_add_comment(of,
                         "nstat (stat:sum): number of beam primaries this run represents; divide tallies derived "
                         "from these particles by it to get a per-primary result");
    mcpl_hdr_add_stat_sum(of, "nstat", (double) nstat);

    memset(&particle, 0, sizeof(particle));
    for (i = 0; i < count; ++i) {
        struct osh_scoring_mcpl_record const *rec = &acc->mcpl_records[i];

        particle.position[0] = rec->position[0];
        particle.position[1] = rec->position[1];
        particle.position[2] = rec->position[2];
        particle.direction[0] = rec->direction[0];
        particle.direction[1] = rec->direction[1];
        particle.direction[2] = rec->direction[2];
        normalize_direction(particle.direction);
        particle.ekin = (rec->ekin > 0.0) ? rec->ekin : 0.0;
        particle.weight = rec->weight;
        particle.pdgcode = rec->pdgcode;
        particle.userflags = rec->gen;
        mcpl_add_particle(of, &particle);
    }

    mcpl_close_outfile(of);
    return OSH_OK;
}
