/**
 * @file osh_transport_prepare.c
 *
 * @brief Build the runtime stopping-power and CSDA-range tables.
 *
 * @details
 * Two paths to fill one (material, projectile, energy) cell:
 *
 *   LOADDEDX path  — material has dedx_table_path != NULL.  The source table
 *                    is loaded (53 log-spaced points typically), then each
 *                    runtime grid point is filled via log-log interpolation.
 *                    Range is then integrated numerically from the resampled SP.
 *
 *   Bethe path     — no LOADDEDX table.  Stopping powers are evaluated directly
 *                    via osh_physics_bethe_eval() at each of the 500 grid points.
 *                    Range integrated the same way.
 *
 * All materials share one energy grid, one projectile list, and one set of grid
 * parameters (log_emin, inv_dlog), so material indices are usable directly as
 * array offsets throughout transport.
 *
 * Range integration:
 *   R(E_n/A) = integral_{E_0/A}^{E_n/A} A * d(E/A) / SP(E/A)
 *   Approximated with the trapezoidal rule on the resampled runtime grid:
 *     dR = 0.5 * (1/SP_i + 1/SP_{i+1}) * A * ((E/A)_{i+1} - (E/A)_i)
 *   Accumulated from the minimum energy upward so that range_csda[k] is the
 *   CSDA range from emin to energy[k].  Runtime transport takes the difference
 *   R(E_in) - R(E_out) for each step.
 */

#include "transport/prepare/osh_transport_prepare.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_const.h"
#include "common/osh_interpolate.h"
#include "common/osh_logger.h"
#include "material/osh_material.h"
#include "material/osh_material_atomic_data.h"
#include "material/osh_material_loaddedx.h"
#include "particle/osh_isotope_db.h"
#include "particle/osh_isotope_db_generated.h"
#include "particle/osh_particle.h"
#include "physics/osh_physics_bethe.h"

/* ---- Helpers -------------------------------------------------------------- */

/**
 * @brief Compute the flat index into the [nmaterials][nprojectiles][nenergy] table.
 */
static inline size_t rt_index(struct osh_transport_runtime_tables const *t, size_t mat, size_t proj, size_t e) {
    return (mat * t->nprojectiles + proj) * t->nenergy + e;
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
static void integrate_range(float const *sp, float *range, double a_proj, struct osh_transport_runtime_tables const *t) {
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
 * @brief Return the default projectile mass number for a given atomic number.
 *
 * @details
 * Reuse the default-isotope selection from the particle database so LOADDEDX
 * source-table import and Bethe fallback agree on which representative ion a
 * projectile column corresponds to.
 */
static unsigned int default_isotope_a(unsigned int z) {
    unsigned int idx_default;

    if (z >= OSH_ISOTOPE_DB_NELEM) {
        return 0u;
    }

    idx_default = osh_isotopes_idx_default[z];
    if (idx_default == OSH_ISOTOPE_DB_ERR || idx_default >= OSH_ISOTOPE_DB_NISO) {
        return 0u;
    }

    return osh_isotope_db[idx_default].a;
}

/**
 * @brief Get the atomic mass to use for one material element.
 *
 * @details
 * Natural elements use the tabulated natural atomic weight from the material
 * atomic-data table. Explicit isotopes use the isotope database. This mirrors
 * material assembly and avoids silently treating explicit isotopes as natural
 * elements during Bethe target preparation.
 */
static double material_element_mass_da(struct material_element const *el) {
    struct isotope iso;
    double mass_da;

    if (!el) {
        return 0.0;
    }

    if (el->a > 0u) {
        if (osh_isotope_from_za(&iso, el->z, el->a)) {
            return iso.amass;
        }
    }

    if (osh_material_natural_atomic_mass_da(el->z, &mass_da) == OSH_OK) {
        return mass_da;
    }

    return 2.0 * (double) el->z;
}

/**
 * @brief Return the nuclear rest mass [MeV/c²] for a projectile (Z, A).
 *
 * @details
 * Looks up the atomic mass from the isotope database and subtracts Z electron
 * masses to obtain the fully-stripped nuclear mass:
 *
 *   M_nuclear = amass [amu] * OSH_AMU  -  Z * m_electron
 *
 * Falls back to A * OSH_AMU (no electron correction) if the isotope is not in
 * the database.  The fallback is accurate to ~0.03% and is only reached for
 * very heavy or exotic nuclei not in the table.
 */
static double projectile_nuclear_mass_mev(unsigned int z, unsigned int a) {
    double mass;

    if (osh_particle_nuclear_mass_mev_from_za(z, a, &mass)) {
        return mass;
    }
    /* Fallback: bare A * u, no electron correction. */
    return (double) a * OSH_AMU;
}

/**
 * @brief Return 1 when every element in a compound has a known MEE value.
 *
 * @details
 * This is the trigger for the SH-compatible element-by-element Bethe path.
 * Compounds with only a material-level mean excitation energy remain on the
 * effective-medium path instead.
 */
static int material_has_complete_element_mee(struct material const *mat) {
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
static void build_bethe_target(struct material const *mat, struct osh_physics_bethe_target *tgt) {
    size_t i;
    double sum_wz_over_a;
    double sum_wz;
    double mee_ln_sum;

    sum_wz_over_a = 0.0;
    sum_wz = 0.0;
    mee_ln_sum = 0.0;

    for (i = 0; i < mat->nelements; ++i) {
        struct material_element const *el;
        double a_i;
        double z_i;
        double w_i;
        double i_i;

        el = &mat->elements[i];
        a_i = material_element_mass_da(el);
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
 * @brief Build one pure-element Bethe target for compound summation.
 *
 * @details
 * This mirrors the old SHIELD-HIT compound treatment more closely than the
 * effective-medium target: each element keeps its own Z, A, density-effect
 * density, and mean excitation energy, and the final mass stopping power is
 * summed by material mass fraction.
 */
static void build_element_bethe_target(struct material_element const *el, struct osh_physics_bethe_target *tgt) {
    struct osh_material_icru_entry entry;
    double mass_da;

    mass_da = material_element_mass_da(el);

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
 */
static void fill_bethe_projectile_column(struct osh_transport_runtime_tables const *t,
                                         size_t mat_idx,
                                         size_t proj_idx,
                                         double dlog,
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
    proj.mass_mev = projectile_nuclear_mass_mev(t->projectile_z[proj_idx], t->projectile_a[proj_idx]);

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
static void fill_bethe_compound_projectile_column(struct osh_transport_runtime_tables const *t,
                                                  size_t mat_idx,
                                                  size_t proj_idx,
                                                  double dlog,
                                                  struct material const *mat) {
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
    proj.mass_mev = projectile_nuclear_mass_mev(t->projectile_z[proj_idx], t->projectile_a[proj_idx]);

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
                                       struct osh_transport_runtime_tables const *t,
                                       size_t mat_idx,
                                       size_t proj_idx) {
    size_t base;
    unsigned int z;
    unsigned int a;
    float sp0;
    float sp1;
    float r1;

    if (osh_log_get_level() > OSH_LOG_DEBUG) {
        return;
    }

    base = rt_index(t, mat_idx, proj_idx, 0);
    z = t->projectile_z[proj_idx];
    a = t->projectile_a[proj_idx];
    sp0 = t->mass_stopping_power[base];
    sp1 = t->mass_stopping_power[base + t->nenergy - 1u];
    r1 = t->range_csda[base + t->nenergy - 1u];

    /*
     * Emit through osh_info() under an explicit debug-level gate. The logger
     * already routes these only for -vv runs via the level check above, and
     * this has proven more reliable in practice for surfacing the runtime
     * table summaries during validate-mode dry runs.
     */
    osh_info("    %s projectile[%zu] Z=%u A=%u: SP(%.3f)=%.6g SP(%.1f)=%.6g Range(%.1f)=%.6g [%s]",
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

enum osh_status osh_transport_prepare(struct material_workspace const *wm,
                                      unsigned int z_max,
                                      struct osh_transport_runtime_tables *tables) {
    size_t nmat, nproj, ne, nbytes_sp, nbytes_range;
    size_t mat_idx, proj_idx, e_idx, base;
    size_t i;
    double dlog, log_e, e_val;
    struct osh_transport_runtime_tables t;
    struct material const *mat;
    struct osh_material_loaddedx_table src;
    struct osh_physics_bethe_target tgt;
    float *sp_col, *rng_col;
    enum osh_status rc;
    if (!wm || !tables)
        return OSH_EINVAL;

    memset(&t, 0, sizeof(t));
    memset(&src, 0, sizeof(src));

    /* ---- Grid parameters ------------------------------------------------- */
    ne = OSH_TRANSPORT_NENERGY;
    t.nenergy = ne;
    t.emin = OSH_TRANSPORT_EMIN;
    t.emax = OSH_TRANSPORT_EMAX;
    t.log_emin = log(OSH_TRANSPORT_EMIN);
    dlog = log(OSH_TRANSPORT_EMAX / OSH_TRANSPORT_EMIN) / (double) (ne - 1u);
    t.inv_dlog = 1.0 / dlog;

    /* ---- Projectile list -------------------------------------------------- */
    /*
     * Expand the caller-requested projectile range to include the widest
     * contiguous LOADDEDX coverage seen in any material. The runtime set is
     * the union of both, so later code can fill higher-Z columns with Bethe
     * fallback when an external table covers fewer projectiles.
     */
    for (i = OSH_MATERIAL_INDEX_FIRST_USER; i < wm->nmaterials; ++i) {
        mat = osh_material_by_index(wm, i);
        if (mat && mat->dedx_table_path) {
            rc = osh_material_loaddedx_table_load(mat->dedx_table_path, &src);
            if (rc != OSH_OK) {
                goto fail;
            }
            if (src.nprojectiles > z_max) {
                z_max = (unsigned int) src.nprojectiles;
            }
            osh_material_loaddedx_table_free(&src);
        }
    }

    nproj = (size_t) z_max;
    if (nproj == 0)
        nproj = 1; /* always at least protons */

    t.nprojectiles = nproj;

    osh_info("Preparing runtime stopping-power tables: %zu energy points from %.3f to %.1f MeV/u", ne, t.emin, t.emax);
    osh_info("Runtime projectile set: %zu representative ions (Z = 1..%zu)", nproj, nproj);

    /* Projectile Z and A arrays. */
    t.projectile_z = calloc(nproj, sizeof(*t.projectile_z));
    t.projectile_a = calloc(nproj, sizeof(*t.projectile_a));
    if (!t.projectile_z || !t.projectile_a) {
        rc = OSH_ENOMEM;
        goto fail;
    }

    for (proj_idx = 0; proj_idx < nproj; ++proj_idx) {
        unsigned int z;
        unsigned int a;

        z = (unsigned int) (proj_idx + 1u);
        a = default_isotope_a(z);
        if (a == 0u) {
            osh_error("Unsupported projectile Z=%u: no default isotope in the isotope database", z);
            rc = OSH_EINVAL;
            goto fail;
        }

        t.projectile_z[proj_idx] = z;
        t.projectile_a[proj_idx] = a;
        osh_info("    projectile[%zu]: Z=%u A=%u", proj_idx, t.projectile_z[proj_idx], t.projectile_a[proj_idx]);
    }

    /* ---- Table storage --------------------------------------------------- */
    nmat = wm->nmaterials; /* includes blackhole (0) and vacuum (1) */
    t.nmaterials = nmat;

    nbytes_sp = nmat * nproj * ne * sizeof(float);
    nbytes_range = nbytes_sp;

    t.mass_stopping_power = calloc(1, nbytes_sp);
    t.range_csda = calloc(1, nbytes_range);
    if (!t.mass_stopping_power || !t.range_csda) {
        rc = OSH_ENOMEM;
        goto fail;
    }

    /* Indices 0 (blackhole) and 1 (vacuum) remain zero. */

    /* ---- Fill per-material rows ------------------------------------------ */
    for (mat_idx = OSH_MATERIAL_INDEX_FIRST_USER; mat_idx < nmat; ++mat_idx) {
        mat = osh_material_by_index(wm, mat_idx);
        if (!mat)
            continue;

        /* ================================================================
         * LOADDEDX path
         * ================================================================ */
        if (mat->dedx_table_path) {
            size_t n_loaddedx;

            rc = osh_material_loaddedx_table_load(mat->dedx_table_path, &src);
            if (rc != OSH_OK)
                goto fail;
            n_loaddedx = src.nprojectiles;

            osh_info("Material '%s': resampling LOADDEDX table %s for %zu projectiles",
                     mat->name,
                     mat->dedx_table_path,
                     n_loaddedx);

            for (proj_idx = 0; proj_idx < nproj && proj_idx < src.nprojectiles; ++proj_idx) {
                base = rt_index(&t, mat_idx, proj_idx, 0);
                sp_col = t.mass_stopping_power + base;

                /* Resample source onto runtime log-uniform grid via log-log. */
                for (e_idx = 0; e_idx < ne; ++e_idx) {
                    log_e = t.log_emin + (double) e_idx * dlog;
                    e_val = exp(log_e);
                    sp_col[e_idx] = (float) osh_interpolate_dloglog(e_val,
                                                                    src.energy_grid,
                                                                    src.mass_stopping_power + proj_idx * src.nenergy,
                                                                    (unsigned int) src.nenergy,
                                                                    OSH_INTERPOLATE_OOB_NEAREST);
                }

                rng_col = t.range_csda + base;
                integrate_range(sp_col, rng_col, (double) t.projectile_a[proj_idx], &t);
                log_runtime_column_summary(mat->name, "LOADDEDX-resampled", &t, mat_idx, proj_idx);
            }

            if (n_loaddedx < nproj) {
                osh_info("Material '%s': filling remaining %zu projectile columns with Bethe fallback",
                         mat->name,
                         nproj - n_loaddedx);
                if (mat->nelements > 1u && material_has_complete_element_mee(mat)) {
                    for (proj_idx = n_loaddedx; proj_idx < nproj; ++proj_idx) {
                        fill_bethe_compound_projectile_column(&t, mat_idx, proj_idx, dlog, mat);
                        log_runtime_column_summary(mat->name, "Bethe-fallback-elements", &t, mat_idx, proj_idx);
                    }
                } else {
                    build_bethe_target(mat, &tgt);
                    for (proj_idx = n_loaddedx; proj_idx < nproj; ++proj_idx) {
                        fill_bethe_projectile_column(&t, mat_idx, proj_idx, dlog, &tgt);
                        log_runtime_column_summary(mat->name, "Bethe-fallback-effective", &t, mat_idx, proj_idx);
                    }
                }
            }

            osh_material_loaddedx_table_free(&src);

            /* ================================================================
             * Bethe path
             * ================================================================ */
        } else {
            if (mat->nelements > 1u && material_has_complete_element_mee(mat)) {
                osh_info("Material '%s': generating Bethe stopping powers for %zu projectiles using element-by-element "
                         "compound mode",
                         mat->name,
                         nproj);
                for (proj_idx = 0; proj_idx < nproj; ++proj_idx) {
                    fill_bethe_compound_projectile_column(&t, mat_idx, proj_idx, dlog, mat);
                    log_runtime_column_summary(mat->name, "Bethe-elements", &t, mat_idx, proj_idx);
                }
            } else {
                osh_info(
                    "Material '%s': generating Bethe stopping powers for %zu projectiles using effective-medium mode",
                    mat->name,
                    nproj);
                build_bethe_target(mat, &tgt);

                for (proj_idx = 0; proj_idx < nproj; ++proj_idx) {
                    fill_bethe_projectile_column(&t, mat_idx, proj_idx, dlog, &tgt);
                    log_runtime_column_summary(mat->name, "Bethe-effective", &t, mat_idx, proj_idx);
                }
            }
        }
    }

    *tables = t;
    return OSH_OK;

fail:
    osh_transport_runtime_tables_free(&t);
    return rc;
}

void osh_transport_runtime_tables_free(struct osh_transport_runtime_tables *tables) {
    if (!tables)
        return;
    free(tables->mass_stopping_power);
    free(tables->range_csda);
    free(tables->projectile_z);
    free(tables->projectile_a);
    memset(tables, 0, sizeof(*tables));
}
