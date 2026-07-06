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

/* Per-step LET context shared by the two-pass LET scorers (DLET/TLET).  Gathered
 * once per step; osh_scoring_estimator_step_let_apply() then resolves the LET a
 * given page scores, honouring that page's Settings medium/density override. */
struct osh_scoring_step_let_ctx {
    int have_proj;       /* projectile found in the SP tables (LET default is table-based) */
    size_t proj_idx;     /* projectile column in the SP tables */
    double e_per_nuc;    /* mean step energy per nucleon [MeV/u] */
    double sp_transport; /* S(transport medium, E) [MeV*cm2/g]; 0 when no table entry */
    double let_default;  /* transport-medium LET [MeV/cm]; de/score_len geometric fallback */
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

/* Gather the transport-medium LET once per step for the two-pass LET scorers.
 * Caller must have rejected neutrals and non-nuclei (z==0 || a==0 || rho<=0) first,
 * so the per-nucleon divide here is safe. */
struct osh_scoring_step_let_ctx osh_scoring_estimator_step_let_gather(struct osh_scoring_runtime const *rt,
                                                                      struct particle const *part,
                                                                      struct step const *st,
                                                                      double score_len);

/* LET [MeV/cm] a page scores for this step: the gathered transport-medium default,
 * or the value under the page's Settings medium/density override when present. */
double osh_scoring_estimator_step_let_apply(struct osh_scoring_step_let_ctx const *ctx,
                                            struct osh_scoring_runtime const *rt,
                                            struct osh_scoring_page_runtime const *page);

/* Flat accumulator index for spatial bin @p spatial_idx at differential offsets
 * (db, db2); inactive differential axes have offset 0.
 *
 * Invariant: page->diff_stride is also the spatial-bin count, because the
 * differential axes are the *outer* strides (a differential axis multiplies the
 * per-spatial-bin count).  The handlers rely on this both here and in their
 * `spatial_idx < page->diff_stride` bound checks; a layout that made a diff axis
 * inner to the spatial bins would break both and must revisit this helper. */
static inline size_t
osh_scoring_estimator_flat_bin(struct osh_scoring_page_runtime const *page, size_t spatial_idx, size_t db, size_t db2) {
    return spatial_idx + (db * page->diff_stride) + (db2 * page->diff2_stride);
}

#endif /* OSH_SCORING_ESTIMATOR_COMMON_H */
