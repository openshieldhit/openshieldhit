#include "physics/neutron/osh_neutron_xsec.h"

#include <math.h>
#include <string.h>

#include "common/osh_diag.h"
#include "openshieldhit/const.h"
#include "physics/neutron/osh_neutron_xsec_data.h"
#include "physics/nuclear/osh_nuclear_tripathi.h"

/*
 * Nuclear radius for the geometric elastic cross-section fallback [cm].
 * σ_geo = π r₀² A^(2/3).  Converted to mb in sigma_geo_mb().
 * r₀ = 1.2 fm matches the standard nuclear-size convention (distinct from
 * Tripathi's r₀ = 1.1 fm, which calibrates the reaction cross section).
 */
#define OSH_NEUTRON_GEO_R0_CM 1.2e-13

/* --------------------------------------------------------------------------
 * Static JEFF-4.0 table
 * -------------------------------------------------------------------------- */

struct xsec_entry {
    int z;
    int a;
    float const *tot_mb;
    float const *el_mb;
    float const *nn_mb;
    float const *n2n_mb;
    float const *ng_mb;
    float const *np_mb;
    float const *na_mb;
};

/*
 * Expand one row of the table from the generated header identifiers.
 * The tag (e.g. "h1", "o16") must match the lowercase condensed nuclide tag
 * produced by tools/condense_neutron_xsec.py.
 */
#define JEFF_ROW(tag, ZZ, AA)                                                                                          \
    {(ZZ),                                                                                                             \
     (AA),                                                                                                             \
     osh_neutron_xsec_##tag##_tot_mb,                                                                                  \
     osh_neutron_xsec_##tag##_el_mb,                                                                                   \
     osh_neutron_xsec_##tag##_nn_mb,                                                                                   \
     osh_neutron_xsec_##tag##_n2n_mb,                                                                                  \
     osh_neutron_xsec_##tag##_ng_mb,                                                                                   \
     osh_neutron_xsec_##tag##_np_mb,                                                                                   \
     osh_neutron_xsec_##tag##_na_mb}

static struct xsec_entry const k_jeff[] = {
    JEFF_ROW(h1, 1, 1),       JEFF_ROW(h2, 1, 2),       JEFF_ROW(he3, 2, 3),      JEFF_ROW(li6, 3, 6),
    JEFF_ROW(li7, 3, 7),      JEFF_ROW(be9, 4, 9),      JEFF_ROW(b10, 5, 10),     JEFF_ROW(b11, 5, 11),
    JEFF_ROW(c12, 6, 12),     JEFF_ROW(n14, 7, 14),     JEFF_ROW(o16, 8, 16),     JEFF_ROW(f19, 9, 19),
    JEFF_ROW(na23, 11, 23),   JEFF_ROW(mg24, 12, 24),   JEFF_ROW(al27, 13, 27),   JEFF_ROW(si28, 14, 28),
    JEFF_ROW(p31, 15, 31),    JEFF_ROW(s32, 16, 32),    JEFF_ROW(cl35, 17, 35),   JEFF_ROW(ar40, 18, 40),
    JEFF_ROW(k39, 19, 39),    JEFF_ROW(ca40, 20, 40),   JEFF_ROW(fe56, 26, 56),   JEFF_ROW(cu63, 29, 63),
    JEFF_ROW(zn64, 30, 64),   JEFF_ROW(zn66, 30, 66),   JEFF_ROW(zn68, 30, 68),   JEFF_ROW(cd113, 48, 113),
    JEFF_ROW(cd114, 48, 114), JEFF_ROW(w182, 74, 182),  JEFF_ROW(w183, 74, 183),  JEFF_ROW(w184, 74, 184),
    JEFF_ROW(w186, 74, 186),  JEFF_ROW(au197, 79, 197), JEFF_ROW(pb208, 82, 208),
};

#undef JEFF_ROW

#define K_JEFF_N ((int) (sizeof(k_jeff) / sizeof(k_jeff[0])))

/* --------------------------------------------------------------------------
 * Interpolation helpers
 * -------------------------------------------------------------------------- */

/*
 * Log-log interpolation of a single cross-section array at energy e_mev.
 * Uses the shared energy grid from the generated header.
 * Linear interpolation is used at threshold crossings (one endpoint zero).
 * Clamps to the edge values outside the grid.
 */
static double interp_one(float const *arr, double e_mev) {
    double le;
    double t;
    double le_lo, le_hi;
    int lo, hi, mid;

    /* clamp at grid edges */
    if (e_mev <= (double) osh_neutron_xsec_egrid_mev[0])
        return (double) arr[0];
    if (e_mev >= (double) osh_neutron_xsec_egrid_mev[OSH_NEUTRON_XSEC_NPOINTS - 1])
        return (double) arr[OSH_NEUTRON_XSEC_NPOINTS - 1];

    /* binary search for bracket [lo, hi] */
    le = log(e_mev);
    lo = 0;
    hi = OSH_NEUTRON_XSEC_NPOINTS - 1;
    while (hi - lo > 1) {
        mid = (lo + hi) / 2;
        if (log((double) osh_neutron_xsec_egrid_mev[mid]) <= le)
            lo = mid;
        else
            hi = mid;
    }

    /* zero-endpoint handling: linear interpolation across threshold */
    if (arr[lo] <= 0.0f || arr[hi] <= 0.0f) {
        if (arr[lo] <= 0.0f && arr[hi] <= 0.0f)
            return 0.0;
        le_lo = (double) osh_neutron_xsec_egrid_mev[lo];
        le_hi = (double) osh_neutron_xsec_egrid_mev[hi];
        t = (e_mev - le_lo) / (le_hi - le_lo);
        return (double) arr[lo] + t * ((double) arr[hi] - (double) arr[lo]);
    }

    /* log-log interpolation */
    le_lo = log((double) osh_neutron_xsec_egrid_mev[lo]);
    le_hi = log((double) osh_neutron_xsec_egrid_mev[hi]);
    t = (le - le_lo) / (le_hi - le_lo);
    return exp(log((double) arr[lo]) + t * (log((double) arr[hi]) - log((double) arr[lo])));
}

/* --------------------------------------------------------------------------
 * Tier-2 optical fallback
 * -------------------------------------------------------------------------- */

static double sigma_geo_mb(int a) {
    double cbrt_a;
    cbrt_a = cbrt((double) a);
    return OSH_M_PI * OSH_NEUTRON_GEO_R0_CM * OSH_NEUTRON_GEO_R0_CM * cbrt_a * cbrt_a * OSH_CM2_TO_MB;
}

static void
lookup_tier2(struct osh_neutron_xsec *xsec, int z, int a, double e_mev, struct osh_neutron_xsec_result *out) {
    double sig_r_cm2;
    double sig_r_mb;
    double sig_el_mb;
    int i;

    /* emit one-time diagnostic per (Z,A) */
    if (xsec->n_warned < OSH_NEUTRON_XSEC_MAX_WARNED) {
        int already;
        already = 0;
        for (i = 0; i < xsec->n_warned; ++i) {
            if (xsec->warned_z[i] == z && xsec->warned_a[i] == a) {
                already = 1;
                break;
            }
        }
        if (!already) {
            OSH_DIAG_INFOF(xsec->diag,
                           "neutron xsec: no JEFF-4.0 data for (Z=%d A=%d), "
                           "using Tripathi + geometric fallback",
                           z,
                           a);
            xsec->warned_z[xsec->n_warned] = z;
            xsec->warned_a[xsec->n_warned] = a;
            xsec->n_warned++;
        }
    }

    /* σ_R [cm²] from Tripathi; neutron projectile: Z=0, A=1 */
    sig_r_cm2 = osh_nuclear_tripathi_sigma(0u, 1u, (double) z, (double) a, e_mev);
    sig_r_mb = sig_r_cm2 * OSH_CM2_TO_MB;

    /* geometric elastic estimate */
    sig_el_mb = sigma_geo_mb(a);

    out->el = sig_el_mb;
    out->tot = sig_r_mb + sig_el_mb;
    out->nonel = sig_r_mb;
    /* sub-channels unknown for Tier-2 — reaction layer treats as generic compound */
    out->nn = 0.0;
    out->n2n = 0.0;
    out->ng = 0.0;
    out->np = 0.0;
    out->na = 0.0;
}

/* --------------------------------------------------------------------------
 * Public API
 * -------------------------------------------------------------------------- */

enum osh_status osh_neutron_xsec_compile(struct osh_diag_sink const *diag, struct osh_neutron_xsec *xsec) {
    if (!xsec)
        return OSH_EINVAL;
    xsec->diag = diag;
    xsec->n_warned = 0;
    return OSH_OK;
}

void osh_neutron_xsec_free(struct osh_neutron_xsec *xsec) {
    /* all data is static; nothing to release */
    (void) xsec;
}

void osh_neutron_xsec_lookup(
    struct osh_neutron_xsec *xsec, int z, int a, double e_mev, struct osh_neutron_xsec_result *out) {
    struct xsec_entry const *entry;
    int i;

    entry = NULL;
    for (i = 0; i < K_JEFF_N; ++i) {
        if (k_jeff[i].z == z && k_jeff[i].a == a) {
            entry = &k_jeff[i];
            break;
        }
    }

    if (!entry) {
        /* -- Tier-2: Tripathi + geometric fallback ----------------------- */
        lookup_tier2(xsec, z, a, e_mev, out);
        return;
    }

    /* -- Tier-1: JEFF-4.0 log-log interpolation ------------------------- */
    out->tot = interp_one(entry->tot_mb, e_mev);
    out->el = interp_one(entry->el_mb, e_mev);
    out->nn = interp_one(entry->nn_mb, e_mev);
    out->n2n = interp_one(entry->n2n_mb, e_mev);
    out->ng = interp_one(entry->ng_mb, e_mev);
    out->np = interp_one(entry->np_mb, e_mev);
    out->na = interp_one(entry->na_mb, e_mev);
    out->nonel = out->tot - out->el;
    if (out->nonel < 0.0)
        out->nonel = 0.0;
}
