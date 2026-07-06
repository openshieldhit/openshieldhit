#ifndef OSH_SCORING_KERNELS_H
#define OSH_SCORING_KERNELS_H

/**
 * @file osh_scoring_kernels.h
 * @brief Pure per-deposit scoring value kernels — the innermost quantity math.
 *
 * @details
 * These are the tiny, branch-free, pointer-free value functions each estimator
 * books.  Everything they need arrives as plain scalars (energy, path length,
 * density, a pre-gathered scaling factor), so they carry no `struct step` / page /
 * geometry state.  That makes them:
 *   - **readable** — the quantity is one line;
 *   - **SoA-ready** — a future batched wavefront can call them per lane / vectorise
 *     them; only the surrounding driver (gather + filter + scatter) changes;
 *   - **unit-testable** — see tests/unit/test_osh_scoring_kernels.c.
 *
 * **Step vs point are deliberately separate.** A step deposit spreads a quantity
 * over a track segment (per-crossing, path-weighted); a point deposit books the
 * whole quantity at one site (no track, no path).  They are not the same function:
 * `step_energy = de·(path/score_len)` but `point_energy = de`; likewise dose.
 * Fluence and LET/Qeff have no point form (they need a track length).
 *
 * The exact float operation order below is chosen to match the historical inline
 * expressions bit-for-bit (e.g. `de·(path/score_len)`, not `(de·path)/score_len`),
 * so this extraction is a pure refactor.
 */

/* ---- Step (track-segment) deposit kernels — booked per crossing ---------- */

/** Energy deposited in a bin along a step [MeV]: fraction path_len/score_len of de. */
static inline double osh_kernel_step_energy(double de, double path_len, double score_len) {
    return de * (path_len / score_len);
}

/** Fluence contribution along a step [cm]: the track length (postprocess ÷volume). */
static inline double osh_kernel_step_fluence(double path_len) {
    return path_len;
}

/** Dose contribution along a step [MeV·cm³/g]: path length × the pre-gathered
 *  `dose_scale = de/(score_len·rho)·SPR` (postprocess ÷volume → MeV/g). */
static inline double osh_kernel_step_dose(double path_len, double dose_scale) {
    return path_len * dose_scale;
}

/** Dose weight of a crossing for the two-pass DLET/DQEFF averages [MeV]. */
static inline double osh_kernel_dose_weight(double de, double path_len, double score_len) {
    return de * path_len / score_len;
}

/** Track weight of a crossing for the two-pass TLET/TQEFF averages [cm]. */
static inline double osh_kernel_track_weight(double ds, double path_len, double score_len) {
    return ds * path_len / score_len;
}

/* ---- Point (single-site) deposit kernels — the whole quantity at one bin -- */

/** Energy of a point deposit [MeV]: the full release, no path fraction. */
static inline double osh_kernel_point_energy(double de) {
    return de;
}

/** Dose of a point deposit [MeV·cm³/g]: energy per mass at the site, with the
 *  dose-to-medium stopping-power ratio (@p sp_ratio = 1 when no override).
 *  Postprocess ÷volume yields MeV/g. */
static inline double osh_kernel_point_dose(double de, double rho, double sp_ratio) {
    return (de / rho) * sp_ratio;
}

#endif /* OSH_SCORING_KERNELS_H */
