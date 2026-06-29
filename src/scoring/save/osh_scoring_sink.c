#include "scoring/save/osh_scoring_sink.h"

#include "scoring/save/osh_scoring_save.h"

/* Native file sink: write the selected outputs to disk via the existing save
 * layer.  Adding const back to the void* ctx is well-defined; we never write
 * through it. */
static enum osh_status file_sink_save(void *ctx,
                                      struct osh_scoring_runtime const *rt,
                                      unsigned long long completed_nstat,
                                      size_t const *want,
                                      size_t n_want) {
    struct osh_scoring_file_sink const *fs = (struct osh_scoring_file_sink const *) ctx;

    if (!fs || !fs->ws) {
        return OSH_EINVAL;
    }
    return osh_scoring_save_outputs(fs->ws, rt, completed_nstat, want, n_want);
}

enum osh_status osh_scoring_file_sink_init(struct osh_scoring_file_sink *fs,
                                           struct osh_scoring_workspace const *ws,
                                           struct osh_scoring_sink *out) {
    if (!fs || !ws || !out) {
        return OSH_EINVAL;
    }
    fs->ws = ws;
    out->save = file_sink_save;
    out->ctx = fs;
    return OSH_OK;
}
