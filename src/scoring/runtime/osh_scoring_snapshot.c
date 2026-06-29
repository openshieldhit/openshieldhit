#include "scoring/runtime/osh_scoring_snapshot.h"

enum osh_status osh_scoring_snapshot_save(struct osh_scoring_sink const *sink,
                                          struct osh_scoring_shadow *shadow,
                                          unsigned long long completed_nstat,
                                          size_t const *want,
                                          size_t n_want) {
    enum osh_status rc;

    if (!sink || !sink->save || !shadow) {
        return OSH_EINVAL;
    }

    rc = osh_scoring_shadow_refresh(shadow);
    if (rc != OSH_OK) {
        return rc;
    }

    return sink->save(sink->ctx, osh_scoring_shadow_view(shadow), completed_nstat, want, n_want);
}
