/**
 * @file osh_material_compile.c
 *
 * @brief Build the runtime stopping-power and CSDA-range tables.
 *
 * @details
 * Two paths to fill one (material, projectile, energy) cell:
 *
 *   Override path  — material has one or more dE/dx overrides. Each covered
 *                    projectile curve is resampled onto the runtime grid via
 *                    log-log interpolation. Range is then integrated
 *                    numerically from the resampled stopping power.
 *
 *   Bethe path     — no material-owned override covers the projectile. Stopping
 *                    powers are evaluated directly via
 *                    osh_physics_bethe_eval() at each of the 500 grid points.
 *                    Range integrated the same way.
 *
 * All materials share one energy grid, one projectile list, and one set of
 * grid parameters (log_emin, inv_dlog), so material indices are usable
 * directly as array offsets throughout transport.
 *
 * Range integration:
 *   R(E_n/A) = integral_{E_0/A}^{E_n/A} A * d(E/A) / SP(E/A)
 *   Approximated with the trapezoidal rule on the resampled runtime grid:
 *     dR = 0.5 * (1/SP_i + 1/SP_{i+1}) * A * ((E/A)_{i+1} - (E/A)_i)
 *   Accumulated from the minimum energy upward so that range_csda[k] is the
 *   CSDA range from emin to energy[k].  Runtime transport takes the difference
 *   R(E_in) - R(E_out) for each step.
 */

#include "material/runtime/osh_material_compile.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_diag.h"
#include "common/osh_interpolate.h"
#include "material/osh_material.h"
#include "material/osh_material_atomic_data.h"
#include "material/osh_material_icru.h"
#include "particle/osh_particle.h"
#include "physics/osh_physics_bethe.h"
#include "voxel/osh_voxel_hu_lut.h"

/* ---- Helpers -------------------------------------------------------------- */

/**
 * @brief Compute the flat index into the [nmaterials][nprojectiles][nenergy] table.
 */
static inline size_t rt_index(struct osh_material_runtime const *t, size_t mat, size_t proj, size_t e) {
    return (mat * t->nprojectiles + proj) * t->nenergy + e;
}

/**
 * @brief Find one exact material-owned dE/dx override by projectile Z.
 */
static struct osh_material_dedx_override const *find_material_dedx_override(struct osh_material const *mat,
                                                                            unsigned int projectile_z) {
    size_t i;

    if (!mat) {
        return NULL;
    }

    i = 0u;
    while (i < mat->ndedx_overrides) {
        if (mat->dedx_overrides[i].projectile_z == projectile_z) {
            return &mat->dedx_overrides[i];
        }
        i++;
    }

    return NULL;
}

/**
 * @brief Log-log interpolation variant for double-valued source curves.
 *
 * @details
 * The common interpolation helper stores legacy table ordinates as float.
 * Material-owned public overrides keep source values in double, so this local
 * helper mirrors the same interpolation rule without widening common API
 * surface in this branch.
 */
static double interpolate_override_dloglog(double xin, double const *xx, double const *ff, size_t n, int mode) {
    long int i;
    double x0, x1, y0, y1, log_x, log_x0, log_x1, t;

    if (xin < xx[0]) {
        switch (mode) {
        case OSH_INTERPOLATE_OOB_ERR:
            return NAN;
        case OSH_INTERPOLATE_OOB_ZERO:
            return 0.0;
        case OSH_INTERPOLATE_OOB_NEAREST:
            return ff[0];
        default:
            x0 = xx[0];
            x1 = xx[1];
            y0 = ff[0];
            y1 = ff[1];
            break;
        }
    } else if (xin > xx[n - 1u]) {
        switch (mode) {
        case OSH_INTERPOLATE_OOB_ERR:
            return NAN;
        case OSH_INTERPOLATE_OOB_ZERO:
            return 0.0;
        case OSH_INTERPOLATE_OOB_NEAREST:
            return ff[n - 1u];
        default:
            x0 = xx[n - 2u];
            x1 = xx[n - 1u];
            y0 = ff[n - 2u];
            y1 = ff[n - 1u];
            break;
        }
    } else {
        i = osh_binary_search_d(xin, xx, (unsigned long int) n);
        if (i < 0) {
            return NAN;
        }
        x0 = xx[i];
        x1 = xx[i + 1];
        y0 = ff[i];
        y1 = ff[i + 1];
    }

    if (x0 <= 0.0 || x1 <= 0.0 || y0 <= 0.0 || y1 <= 0.0) {
        t = (xin - x0) / (x1 - x0);
        return y0 + t * (y1 - y0);
    }

    log_x = log(xin);
    log_x0 = log(x0);
    log_x1 = log(x1);
    t = (log_x - log_x0) / (log_x1 - log_x0);
    return exp(log(y0) + t * (log(y1) - log(y0)));
}

/**
 * @brief Integrate CSDA range from stopping powers on the runtime grid.
 *
 * @details
 * Trapezoidal rule on the (E, 1/SP) curve, accumulated from k=0 upward:
 *
 *   range[k] = range[k-1] + 0.5*(1/sp[k-1] + 1/sp[k]) * (E[k] - E[k-1])
 *
 * Energy grid is not stored explicitly; each point is exp(log_emin + k*dlog).
 * We recompute the two adjacent energies per step to avoid a 500-element
 * auxiliary array.
 *
 * @param[in]  sp       Stopping power array [nenergy] [MeV cm²/g].
 * @param[out] range    Range array [nenergy] [g/cm²], range[0] = 0.
 * @param[in]  a_proj   Projectile mass number A.
 * @param[in]  t        Runtime tables (grid parameters, nenergy).
 */
static void integrate_range(float const *sp, float *range, double a_proj, struct osh_material_runtime const *t) {
    size_t k;
    double dlog = 1.0 / t->inv_dlog;
    double log_e_prev, log_e_curr, e_prev, e_curr, inv_sp_prev, inv_sp_curr, r;

    range[0] = 0.0f;
    log_e_prev = t->log_emin;
    e_prev = t->emin;
    inv_sp_prev = (sp[0] > 0.0f) ? (1.0 / (double) sp[0]) : 0.0;
    r = 0.0;

    for (k = 1; k < t->nenergy; ++k) {
        log_e_curr = t->log_emin + (double) k * dlog;
        e_curr = exp(log_e_curr);
        inv_sp_curr = (sp[k] > 0.0f) ? (1.0 / (double) sp[k]) : 0.0;

        r += 0.5 * (inv_sp_prev + inv_sp_curr) * a_proj * (e_curr - e_prev);
        range[k] = (float) r;

        log_e_prev = log_e_curr;
        e_prev = e_curr;
        inv_sp_prev = inv_sp_curr;
    }

    (void) log_e_prev; /* suppress unused-variable warning */
}

/* ---- Projectile and Bethe helpers ----------------------------------------- */

/**
 * @brief Return 1 when every element in a compound has a known MEE value.
 *
 * @details
 * This is the trigger for the SH-compatible element-by-element Bethe path.
 * Compounds with only a material-level mean excitation energy remain on the
 * effective-medium path instead.
 */
static int material_has_complete_element_mee(struct osh_material const *mat) {
    size_t i;

    if (!mat || mat->nelements == 0u) {
        return 0;
    }

    for (i = 0; i < mat->nelements; ++i) {
        if (mat->elements[i].mean_excitation_energy <= 0.0) {
            return 0;
        }
    }

    return 1;
}

/**
 * @brief Build the effective Bethe target descriptor for one material.
 *
 * @details
 * The Bethe model consumes the material only through Z/A, density, and the
 * material-level mean excitation energy. For compounds with explicit isotope
 * components we therefore use the correct isotope mass in Z/A, but still pass
 * one effective target descriptor into the pure Bethe evaluator.
 */
static void build_bethe_target(struct osh_material const *mat, struct osh_physics_bethe_target *tgt) {
    size_t i;
    double sum_wz_over_a;
    double sum_wz;
    double mee_ln_sum;

    sum_wz_over_a = 0.0;
    sum_wz = 0.0;
    mee_ln_sum = 0.0;

    for (i = 0; i < mat->nelements; ++i) {
        struct osh_material_element const *el;
        double a_i;
        double z_i;
        double w_i;
        double i_i;

        el = &mat->elements[i];
        if (osh_material_atomic_mass_da(el->z, el->a, &a_i) != OSH_OK) {
            a_i = 2.0 * (double) el->z;
        }
        z_i = (double) el->z;
        w_i = (el->mass_fraction >= 0.0) ? el->mass_fraction : 0.0;

        if (a_i <= 0.0) {
            continue;
        }

        sum_wz_over_a += w_i * z_i / a_i;
        sum_wz += w_i * z_i;

        if (el->mean_excitation_energy > 0.0) {
            i_i = el->mean_excitation_energy;
            mee_ln_sum += w_i * (z_i / a_i) * log(i_i);
        }
    }

    if (sum_wz_over_a <= 0.0) {
        sum_wz_over_a = 1.0;
    }
    if (sum_wz <= 0.0) {
        sum_wz = 1.0;
    }

    tgt->z_mean = sum_wz;
    tgt->a_mean = sum_wz / sum_wz_over_a;
    tgt->rho = (mat->rho > 0.0) ? mat->rho : 1.0;

    /*
     * Material-level mean excitation energy is authoritative for transport.
     * Only if it is absent do we derive a Bragg-additive compound value from
     * explicit element-level entries.
     */
    if (mat->mean_excitation_energy > 0.0) {
        tgt->i_value = mat->mean_excitation_energy;
    } else if (mee_ln_sum != 0.0) {
        tgt->i_value = exp(mee_ln_sum / sum_wz_over_a);
    } else {
        tgt->i_value = 78.0;
    }
}

/**
 * @brief Compute per-material atomic scalars used by MCS and straggling.
 *
 * @details
 * Fills three scalars derived from the element composition:
 *
 *   z_mean   = sum_i(w_i * Z_i)          effective atomic number
 *   z_over_a = sum_i(w_i * Z_i / A_i)    Z/A [mol/g] for Bohr straggling
 *   x0_gcm2  = 1 / sum_i(w_i / X0_i)     radiation length [g/cm²]
 *              where X0_i = 716.4 * A_i / [Z_i*(Z_i+1)*ln(287/√Z_i)]
 *              is the approximate element radiation length (PDG eq. 34.25).
 *
 * The radiation-length formula is the PDG first-order approximation and gives
 * ~2% accuracy for most elements.  Compounds use the standard mixture rule.
 * Returns 0 for all three on invalid input (zero atoms or missing mass data).
 */
static void
compute_material_atomic(struct osh_material const *mat, float *z_mean_out, float *z_over_a_out, float *x0_gcm2_out) {
    size_t i;
    double sum_wz;
    double sum_wz_over_a;
    double inv_x0;
    double z_i;
    double a_i;
    double w_i;
    double x0_i;

    sum_wz = 0.0;
    sum_wz_over_a = 0.0;
    inv_x0 = 0.0;

    if (!mat || mat->nelements == 0u) {
        *z_mean_out = 0.0f;
        *z_over_a_out = 0.0f;
        *x0_gcm2_out = 0.0f;
        return;
    }

    for (i = 0; i < mat->nelements; ++i) {
        z_i = (double) mat->elements[i].z;
        w_i = (mat->elements[i].mass_fraction >= 0.0) ? mat->elements[i].mass_fraction : 0.0;

        if (osh_material_atomic_mass_da(mat->elements[i].z, mat->elements[i].a, &a_i) != OSH_OK) {
            a_i = 2.0 * z_i;
        }
        if (a_i <= 0.0 || z_i < 1.0) {
            continue;
        }

        sum_wz += w_i * z_i;
        sum_wz_over_a += w_i * z_i / a_i;

        /* PDG approximate radiation length for pure element [g/cm²]
         * X0 = 716.4 * A / [Z*(Z+1)*ln(287/√Z)]  (PDG eq. 34.25) */
        x0_i = 716.4 * a_i / (z_i * (z_i + 1.0) * log(287.0 / sqrt(z_i)));
        if (x0_i > 0.0) {
            inv_x0 += w_i / x0_i;
        }
    }

    *z_mean_out = (float) sum_wz;
    *z_over_a_out = (float) sum_wz_over_a;
    *x0_gcm2_out = (inv_x0 > 0.0) ? (float) (1.0 / inv_x0) : 0.0f;
}

/**
 * @brief Build one pure-element Bethe target for compound summation.
 *
 * @details
 * This mirrors the old SHIELD-HIT compound treatment more closely than the
 * effective-medium target: each element keeps its own Z, A, density-effect
 * density, and mean excitation energy, and the final mass stopping power is
 * summed by material mass fraction.
 */
static void build_element_bethe_target(struct osh_material_element const *el, struct osh_physics_bethe_target *tgt) {
    struct osh_material_icru_entry entry;
    double mass_da;

    if (osh_material_atomic_mass_da(el->z, el->a, &mass_da) != OSH_OK) {
        mass_da = 2.0 * (double) el->z;
    }

    tgt->z_mean = (double) el->z;
    tgt->a_mean = (mass_da > 0.0) ? mass_da : 2.0 * (double) el->z;
    tgt->i_value = el->mean_excitation_energy;

    if (osh_material_icru_lookup((int) el->z, &entry) == OSH_OK && entry.rho > 0.0) {
        tgt->rho = entry.rho;
    } else {
        tgt->rho = 1.0;
    }
}

/**
 * @brief Fill one runtime projectile column with Bethe stopping power and range.
 *
 * @returns OSH_OK on success, or OSH_EINVAL if projectile mass data is missing.
 */
static enum osh_status fill_bethe_projectile_column(struct osh_material_runtime const *t,
                                                    size_t mat_idx,
                                                    size_t proj_idx,
                                                    double dlog,
                                                    struct osh_diag_sink const *diag,
                                                    struct osh_physics_bethe_target const *tgt) {
    struct osh_physics_bethe_projectile proj;
    struct osh_physics_bethe_sewn sewn;
    float *sp_col;
    float *rng_col;
    size_t base;
    size_t e_idx;
    double log_e;
    double e_val;
    double sp_val;

    proj.z = (double) t->projectile_z[proj_idx];
    proj.a = (double) t->projectile_a[proj_idx];
    if (!osh_particle_nuclear_mass_mev_from_za(t->projectile_z[proj_idx], t->projectile_a[proj_idx], &proj.mass_mev)) {
        OSH_DIAG_ERRORF(diag,
                        "transport: missing nuclear rest mass for projectile Z=%u A=%u",
                        t->projectile_z[proj_idx],
                        t->projectile_a[proj_idx]);
        return OSH_EINVAL;
    }

    osh_physics_bethe_sewn_compute(&proj, tgt, &sewn);

    base = rt_index(t, mat_idx, proj_idx, 0);
    sp_col = t->mass_stopping_power + base;

    for (e_idx = 0; e_idx < t->nenergy; ++e_idx) {
        log_e = t->log_emin + (double) e_idx * dlog;
        e_val = exp(log_e);
        sp_val = osh_physics_bethe_eval(e_val, &proj, tgt, &sewn);
        sp_col[e_idx] = (float) (sp_val > 0.0 ? sp_val : 0.0);
    }

    rng_col = t->range_csda + base;
    integrate_range(sp_col, rng_col, proj.a, t);
    return OSH_OK;
}

/**
 * @brief Fill one runtime projectile column by summing elemental Bethe terms.
 *
 * @details
 * For compounds with element-level MEE values available, this follows the old
 * SHIELD-HIT style more closely than the effective-medium approximation:
 *
 *   S_material = sum_i w_i * S_i
 *
 * where w_i is the element mass fraction and S_i is the pure-element mass
 * stopping power evaluated with that element's own Z, A, density and I.
 */
static enum osh_status fill_bethe_compound_projectile_column(struct osh_material_runtime const *t,
                                                             size_t mat_idx,
                                                             size_t proj_idx,
                                                             double dlog,
                                                             struct osh_diag_sink const *diag,
                                                             struct osh_material const *mat) {
    struct osh_physics_bethe_projectile proj;
    struct osh_physics_bethe_target *elt_tgts;
    struct osh_physics_bethe_sewn *elt_sewns;
    float *sp_col;
    float *rng_col;
    size_t base;
    size_t e_idx;
    size_t i;
    proj.z = (double) t->projectile_z[proj_idx];
    proj.a = (double) t->projectile_a[proj_idx];
    if (!osh_particle_nuclear_mass_mev_from_za(t->projectile_z[proj_idx], t->projectile_a[proj_idx], &proj.mass_mev)) {
        OSH_DIAG_ERRORF(diag,
                        "transport: missing nuclear rest mass for projectile Z=%u A=%u",
                        t->projectile_z[proj_idx],
                        t->projectile_a[proj_idx]);
        return OSH_EINVAL;
    }

    base = rt_index(t, mat_idx, proj_idx, 0);
    sp_col = t->mass_stopping_power + base;

    elt_tgts = (struct osh_physics_bethe_target *) malloc(mat->nelements * sizeof(*elt_tgts));
    elt_sewns = (struct osh_physics_bethe_sewn *) malloc(mat->nelements * sizeof(*elt_sewns));

    if (elt_tgts && elt_sewns) {
        for (i = 0; i < mat->nelements; ++i) {
            build_element_bethe_target(&mat->elements[i], &elt_tgts[i]);
            osh_physics_bethe_sewn_compute(&proj, &elt_tgts[i], &elt_sewns[i]);
        }
    }

    for (e_idx = 0; e_idx < t->nenergy; ++e_idx) {
        double log_e;
        double e_val;
        double sp_sum;

        log_e = t->log_emin + (double) e_idx * dlog;
        e_val = exp(log_e);
        sp_sum = 0.0;

        for (i = 0; i < mat->nelements; ++i) {
            double sp_i;

            if (elt_tgts && elt_sewns) {
                sp_i = osh_physics_bethe_eval(e_val, &proj, &elt_tgts[i], &elt_sewns[i]);
            } else {
                struct osh_physics_bethe_target elt_tgt;
                struct osh_physics_bethe_sewn elt_sewn;

                /* Fallback keeps correctness if temporary precompute buffers
                 * cannot be allocated during setup.
                 */
                build_element_bethe_target(&mat->elements[i], &elt_tgt);
                osh_physics_bethe_sewn_compute(&proj, &elt_tgt, &elt_sewn);
                sp_i = osh_physics_bethe_eval(e_val, &proj, &elt_tgt, &elt_sewn);
            }
            if (sp_i > 0.0) {
                sp_sum += mat->elements[i].mass_fraction * sp_i;
            }
        }

        sp_col[e_idx] = (float) sp_sum;
    }

    free(elt_tgts);
    free(elt_sewns);

    rng_col = t->range_csda + base;
    integrate_range(sp_col, rng_col, proj.a, t);
    return OSH_OK;
}

/**
 * @brief Print one runtime-column summary at debug level.
 *
 * @details
 * The runtime tables themselves are large and not suitable for full dumps
 * during validation. A compact column summary makes it visible that the
 * resampled/generated stopping-power and range arrays were actually built:
 * first and last SP values, and the final CSDA range at emax.
 */
static void log_runtime_column_summary(char const *material_name,
                                       char const *source_tag,
                                       struct osh_material_runtime const *t,
                                       size_t mat_idx,
                                       size_t proj_idx,
                                       struct osh_diag_sink const *diag) {
    size_t base;
    unsigned int z;
    unsigned int a;
    float sp0;
    float sp1;
    float r1;

    if (!diag || !diag->emit || diag->min_level > OSH_DIAG_LEVEL_DEBUG) {
        return;
    }

    base = rt_index(t, mat_idx, proj_idx, 0);
    z = t->projectile_z[proj_idx];
    a = t->projectile_a[proj_idx];
    sp0 = t->mass_stopping_power[base];
    sp1 = t->mass_stopping_power[base + t->nenergy - 1u];
    r1 = t->range_csda[base + t->nenergy - 1u];

    /*
     * Emit through the diagnostics sink under an explicit debug-level gate.
     * This has proven useful for surfacing runtime-table summaries during
     * validate-mode dry runs without adding noise to the hot path.
     */
    OSH_DIAG_INFOF(diag,
                   "    %s projectile[%zu] Z=%u A=%u: SP(%.3f)=%.6g SP(%.1f)=%.6g Range(%.1f)=%.6g [%s]",
                   material_name ? material_name : "(unnamed)",
                   proj_idx,
                   z,
                   a,
                   t->emin,
                   sp0,
                   t->emax,
                   sp1,
                   t->emax,
                   r1,
                   source_tag ? source_tag : "runtime");
}

/* ---- Public API ------------------------------------------------------------ */

enum osh_status osh_material_compile(struct osh_material_workspace const *wm,
                                     unsigned int z_max,
                                     struct osh_diag_sink const *diag,
                                     struct osh_material_runtime *tables) {
    size_t nmat, nproj, ne, nbytes_sp, nbytes_range;
    size_t mat_idx, proj_idx, e_idx, base;
    size_t i, j;
    unsigned int z, a;
    double dlog, log_e, e_val;
    double mass_mev;
    struct osh_material_runtime t;
    struct osh_material const *mat;
    struct osh_material_dedx_override const *dedx_override;
    struct osh_physics_bethe_target tgt;
    float *sp_col, *rng_col;
    int use_compound_bethe;
    enum osh_status rc;
    if (!wm || !tables)
        return OSH_EINVAL;

    memset(&t, 0, sizeof(t));

    /* ---- Grid parameters ------------------------------------------------- */
    ne = OSH_MATERIAL_RUNTIME_NENERGY;
    t.nenergy = ne;
    t.emin = OSH_MATERIAL_RUNTIME_EMIN;
    t.emax = OSH_MATERIAL_RUNTIME_EMAX;
    t.log_emin = log(OSH_MATERIAL_RUNTIME_EMIN);
    dlog = log(OSH_MATERIAL_RUNTIME_EMAX / OSH_MATERIAL_RUNTIME_EMIN) / (double) (ne - 1u);
    t.inv_dlog = 1.0 / dlog;

    /* ---- Projectile list -------------------------------------------------- */
    /* Expand the caller-requested projectile range to include every projectile
     * Z mentioned by any material dE/dx override. */
    for (i = OSH_MATERIAL_INDEX_FIRST_USER; i < wm->nmaterials; ++i) {
        mat = osh_material_by_index(wm, i);
        if (!mat) {
            continue;
        }
        for (j = 0u; j < mat->ndedx_overrides; ++j) {
            if (mat->dedx_overrides[j].projectile_z > z_max) {
                z_max = mat->dedx_overrides[j].projectile_z;
            }
        }
    }

    nproj = (size_t) z_max;
    if (nproj == 0)
        nproj = 1; /* always at least protons */

    t.nprojectiles = nproj;

    OSH_DIAG_INFOF(
        diag, "Preparing runtime stopping-power tables: %zu energy points from %.3f to %.1f MeV/u", ne, t.emin, t.emax);
    OSH_DIAG_INFOF(diag, "Runtime projectile set: %zu representative ions (Z = 1..%zu)", nproj, nproj);

    /* Projectile Z, A, and nuclear-mass arrays. */
    t.projectile_z = calloc(nproj, sizeof(*t.projectile_z));
    t.projectile_a = calloc(nproj, sizeof(*t.projectile_a));
    t.projectile_mass_mev = calloc(nproj, sizeof(*t.projectile_mass_mev));
    if (!t.projectile_z || !t.projectile_a || !t.projectile_mass_mev) {
        rc = OSH_ENOMEM;
        goto fail;
    }

    for (proj_idx = 0; proj_idx < nproj; ++proj_idx) {
        z = (unsigned int) (proj_idx + 1u);
        if (!osh_particle_default_isotope_a(z, &a)) {
            OSH_DIAG_ERRORF(diag, "Unsupported projectile Z=%u: no default isotope in the isotope database", z);
            rc = OSH_EINVAL;
            goto fail;
        }

        if (!osh_particle_nuclear_mass_mev_from_za(z, a, &mass_mev)) {
            OSH_DIAG_ERRORF(diag, "Missing nuclear rest mass for projectile Z=%u A=%u", z, a);
            rc = OSH_EINVAL;
            goto fail;
        }

        t.projectile_z[proj_idx] = z;
        t.projectile_a[proj_idx] = a;
        t.projectile_mass_mev[proj_idx] = mass_mev;
        OSH_DIAG_INFOF(diag, "    projectile[%zu]: Z=%u A=%u mass=%.4f MeV", proj_idx, z, a, mass_mev);
    }

    /* ---- Table storage --------------------------------------------------- */
    nmat = wm->nmaterials; /* includes blackhole (0) and vacuum (1) */
    t.nmaterials = nmat;

    nbytes_sp = nmat * nproj * ne * sizeof(float);
    nbytes_range = nbytes_sp;

    t.mass_stopping_power = calloc(1, nbytes_sp);
    t.range_csda = calloc(1, nbytes_range);
    t.rho = calloc(nmat, sizeof(*t.rho));
    t.z_mean = calloc(nmat, sizeof(*t.z_mean));
    t.z_over_a = calloc(nmat, sizeof(*t.z_over_a));
    t.rad_length = calloc(nmat, sizeof(*t.rad_length));
    if (!t.mass_stopping_power || !t.range_csda || !t.rho || !t.z_mean || !t.z_over_a || !t.rad_length) {
        rc = OSH_ENOMEM;
        goto fail;
    }

    /* Indices 0 (blackhole) and 1 (vacuum) remain zero for all arrays. */

    /* ---- Fill per-material rows ------------------------------------------ */
    for (mat_idx = OSH_MATERIAL_INDEX_FIRST_USER; mat_idx < nmat; ++mat_idx) {
        mat = osh_material_by_index(wm, mat_idx);
        if (!mat)
            continue;

        if (mat->rho > 0.0) {
            t.rho[mat_idx] = (float) mat->rho;
        }

        /* Compute per-material atomic scalars (z_mean, z/a, X₀) for MCS and
         * straggling.  These depend only on composition, not on projectile or
         * energy, so they are computed once here and stored in the flat arrays. */
        compute_material_atomic(mat, &t.z_mean[mat_idx], &t.z_over_a[mat_idx], &t.rad_length[mat_idx]);
        OSH_DIAG_INFOF(diag,
                       "Material '%s': z_mean=%.2f  z/a=%.5f  X0=%.2f g/cm2",
                       mat->name,
                       (double) t.z_mean[mat_idx],
                       (double) t.z_over_a[mat_idx],
                       (double) t.rad_length[mat_idx]);

        /* ================================================================
         * Override and Bethe paths
         * ================================================================ */
        use_compound_bethe = (mat->nelements > 1u && material_has_complete_element_mee(mat)) ? 1 : 0;
        if (!use_compound_bethe) {
            build_bethe_target(mat, &tgt);
        }
        if (mat->ndedx_overrides > 0u) {
            OSH_DIAG_INFOF(diag,
                           "Material '%s': resampling %zu material-owned dE/dx override curves",
                           mat->name,
                           mat->ndedx_overrides);
        } else if (use_compound_bethe) {
            OSH_DIAG_INFOF(
                diag,
                "Material '%s': generating Bethe stopping powers for %zu projectiles using element-by-element "
                "compound mode",
                mat->name,
                nproj);
        } else {
            OSH_DIAG_INFOF(
                diag,
                "Material '%s': generating Bethe stopping powers for %zu projectiles using effective-medium mode",
                mat->name,
                nproj);
        }

        for (proj_idx = 0; proj_idx < nproj; ++proj_idx) {
            dedx_override = find_material_dedx_override(mat, t.projectile_z[proj_idx]);
            if (dedx_override) {
                base = rt_index(&t, mat_idx, proj_idx, 0);
                sp_col = t.mass_stopping_power + base;

                /* Resample source curve onto runtime log-uniform grid via log-log. */
                for (e_idx = 0; e_idx < ne; ++e_idx) {
                    log_e = t.log_emin + (double) e_idx * dlog;
                    e_val = exp(log_e);
                    sp_col[e_idx] = (float) interpolate_override_dloglog(e_val,
                                                                         dedx_override->energy_mev_per_u,
                                                                         dedx_override->dedx_mev_cm2_per_g,
                                                                         dedx_override->npoints,
                                                                         OSH_INTERPOLATE_OOB_NEAREST);
                }

                rng_col = t.range_csda + base;
                integrate_range(sp_col, rng_col, (double) t.projectile_a[proj_idx], &t);
                log_runtime_column_summary(mat->name, "override-resampled", &t, mat_idx, proj_idx, diag);
            } else {
                if (use_compound_bethe) {
                    rc = fill_bethe_compound_projectile_column(&t, mat_idx, proj_idx, dlog, diag, mat);
                    if (rc != OSH_OK) {
                        goto fail;
                    }
                    log_runtime_column_summary(mat->name, "Bethe-elements", &t, mat_idx, proj_idx, diag);
                } else {
                    rc = fill_bethe_projectile_column(&t, mat_idx, proj_idx, dlog, diag, &tgt);
                    if (rc != OSH_OK) {
                        goto fail;
                    }
                    log_runtime_column_summary(mat->name, "Bethe-effective", &t, mat_idx, proj_idx, diag);
                }
            }
        }
    }

    switch (wm->hu_table_type) {
    case OSH_HU_TABLE_NONE:
        break;
    case OSH_HU_TABLE_SCHNEIDER:
    case OSH_HU_TABLE_PERMATASSARI:
        /* Build the HU→density LUT; used by osh_material_runtime_get_rho() for voxel zones. */
        t.hu_rho_lut = (float *) malloc(OSH_VOXEL_HU_LUT_SIZE * sizeof(float));
        if (!t.hu_rho_lut) {
            rc = OSH_ENOMEM;
            goto fail;
        }
        if (wm->hu_table_type == OSH_HU_TABLE_PERMATASSARI) {
            osh_voxel_build_hu_rho_lut_permatassari2020(t.hu_rho_lut);
        } else {
            osh_voxel_build_hu_rho_lut_schneider2000(t.hu_rho_lut);
        }
        break;
    default:
        rc = OSH_EINVAL;
        goto fail;
    }

    *tables = t;
    return OSH_OK;

fail:
    osh_material_runtime_free(&t);
    return rc;
}

void osh_material_runtime_free(struct osh_material_runtime *tables) {
    if (!tables)
        return;
    free(tables->mass_stopping_power);
    free(tables->range_csda);
    free(tables->projectile_z);
    free(tables->projectile_a);
    free(tables->projectile_mass_mev);
    free(tables->rho);
    free(tables->z_mean);
    free(tables->z_over_a);
    free(tables->rad_length);
    free(tables->hu_rho_lut);
    memset(tables, 0, sizeof(*tables));
}
