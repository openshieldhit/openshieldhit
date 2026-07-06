#ifndef OSH_SCORING_ESTIMATOR_COMMON_H
#define OSH_SCORING_ESTIMATOR_COMMON_H

#include <stddef.h>

#include "common/osh_step.h"
#include "common/raytrace/osh_raytrace.h"
#include "material/runtime/osh_material_runtime.h"
#include "particle/osh_particle.h"
#include "scoring/runtime/osh_scoring_accumulator.h"
#include "scoring/runtime/osh_scoring_runtime.h"

/*
 * Shared estimator helpers.
 *
 * These are internal to src/scoring/runtime.  The step and point estimator files
 * use them after geometry has already produced crossings or a located bin.
 */

struct osh_scoring_dose_sp_ctx {
    int have_proj;    /* projectile has an SP-table entry in the transport medium */
    size_t proj_idx;  /* projectile index into the SP tables */
    double sp_tr;     /* S(transport medium, E) [MeV*cm2/g] */
    double e_per_nuc; /* mean step energy per nucleon [MeV/u] */
};

int osh_scoring_estimator_find_proj_idx(struct osh_material_runtime const *tables,
                                        unsigned int z,
                                        size_t *proj_idx_out);

double osh_scoring_estimator_particle_beta(double e_kin_mev, double rest_mass_mev);
double osh_scoring_estimator_particle_qeff(int z, double beta);

int osh_scoring_estimator_resolve_diff_bins(struct osh_scoring_page_runtime const *page,
                                            struct osh_scoring_runtime const *rt,
                                            struct particle const *part,
                                            struct step const *st,
                                            size_t *db,
                                            size_t *db2);

struct osh_scoring_dose_sp_ctx osh_scoring_estimator_dose_sp_gather(struct osh_scoring_runtime const *rt,
                                                                    struct particle const *part,
                                                                    struct step const *st);

double osh_scoring_estimator_dose_sp_ratio(struct osh_scoring_dose_sp_ctx const *c,
                                           struct osh_scoring_runtime const *rt,
                                           struct osh_scoring_page_runtime const *page);

/* Flat accumulator index for spatial bin @p spatial_idx at differential offsets
 * (db, db2); inactive differential axes have offset 0. */
static inline size_t
osh_scoring_estimator_flat_bin(struct osh_scoring_page_runtime const *page, size_t spatial_idx, size_t db, size_t db2) {
    return spatial_idx + (db * page->diff_stride) + (db2 * page->diff2_stride);
}

#endif /* OSH_SCORING_ESTIMATOR_COMMON_H */
