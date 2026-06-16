/*
 * Cylindrical (R,Z) voxel traversal using the alpha-parametric formulation of
 * Jacobs (1998), adapted from the original SHIELD-HIT12A sh_jacobs_cyl.c by
 * Niels Bassler / David C. Hansen.
 *
 * Grid field convention (see osh_raytrace_cyl.h):
 *   grid->origin[0], spacing[0], n[0]  = r_min, dr, nr
 *   grid->origin[2], spacing[2], n[2]  = z_min, dz, nz
 *   index [1] unused.
 *
 * Flat index: ir + nr*iz  (R varies fastest).
 * Max crossings: 2*nr + nz.
 *
 * Reference: F. Jacobs, E. Sundermann, B. De Sutter, M. Christiaens,
 * I. Lemahieu, "A fast algorithm to calculate the exact radiological path
 * through a pixel or voxel space", J Comput Inf Technol. 1998;6(1):89-94.
 * https://hrcak.srce.hr/150245
 */

#include <math.h>

#include "common/raytrace/osh_raytrace_cyl.h"

#define CYL_EPS 1.0e-9
#define CYL_EPSINV 1.0e9

#define CYL_MIN(a, b) (((a) < (b)) ? (a) : (b))
#define CYL_MAX(a, b) (((a) > (b)) ? (a) : (b))
#define CYL_MIN3(a, b, c) (CYL_MIN(CYL_MIN(a, b), c))
#define CYL_MAX3(a, b, c) (CYL_MAX(CYL_MAX(a, b), c))
#define CYL_MIN4(a, b, c, d) (CYL_MIN(CYL_MIN(a, b), CYL_MIN(c, d)))

/* ---- small helpers ------------------------------------------------------- */

static inline int valid_index(size_t ir, size_t iz, size_t nr, size_t nz) {
    return ir < nr && iz < nz;
}

static inline size_t flatten_idx(size_t ir, size_t iz, size_t nr) {
    return ir + nr * iz;
}

/* Determinant for quadratic circle intersection: d = gb²-4·ga·(gc-r²).
 * Positive means two real intersections exist. */
static inline double det2(double ga, double gb, double gc, double r2) {
    return gb * gb - 4.0 * ga * (gc - r2);
}

/* Solve ga·α²+gb·α+(gc-r²)=0; sets a1,a2=CYL_EPSINV if no real roots. */
static inline void slv2(double ga, double gb, double gc, double r2, double *a1, double *a2) {
    double d;  /* discriminant */
    double sq; /* half-spread: sqrt(d)/(2*ga) */
    double bb; /* midpoint: -gb/(2*ga) */

    d = det2(ga, gb, gc, r2);
    if (d <= 0.0) {
        *a1 = CYL_EPSINV;
        *a2 = CYL_EPSINV;
        return;
    }
    sq = sqrt(d) / (2.0 * ga);
    bb = -gb / (2.0 * ga);
    *a1 = bb - sq;
    *a2 = bb + sq;
}

/* Alpha of the next Z-plane boundary after current bin iz. */
static inline double next_alpha_z(double const ds[3], double pz, double z_min, double dz, size_t iz) {
    double z1; /* Z coordinate of lower bin face */
    double z2; /* Z coordinate of upper bin face */
    double a1; /* alpha for z1 crossing */
    double a2; /* alpha for z2 crossing */

    z1 = z_min + dz * (double) iz;
    z2 = z_min + dz * (double) (iz + 1u);
    a1 = (z1 - pz) / ds[2];
    a2 = (z2 - pz) / ds[2];
    return CYL_MAX(a1, a2);
}

/* Alpha of the next R-shell boundary after current bin ir, from current alpha.
 * Checks both the inner AND outer shell of bin ir so a turning-point ray (one
 * whose r decreases to a minimum then increases) is handled correctly: the
 * outward crossing of the outer shell is visible even when the ray is currently
 * moving inward.  No direction heuristic is needed. */
static inline double next_alpha_r(double ga, double gb, double gc, double r_min, double dr, size_t ir, double alpha) {
    double r_inner;  /* radius of inner shell boundary of bin ir */
    double r_outer;  /* radius of outer shell boundary of bin ir */
    double a11, a12; /* roots for r_inner circle */
    double a21, a22; /* roots for r_outer circle */

    r_inner = r_min + dr * (double) ir;
    r_outer = r_min + dr * (double) (ir + 1u);
    slv2(ga, gb, gc, r_inner * r_inner, &a11, &a12);
    slv2(ga, gb, gc, r_outer * r_outer, &a21, &a22);
    if (a11 <= alpha)
        a11 = CYL_EPSINV;
    if (a12 <= alpha)
        a12 = CYL_EPSINV;
    if (a21 <= alpha)
        a21 = CYL_EPSINV;
    if (a22 <= alpha)
        a22 = CYL_EPSINV;
    return CYL_MIN4(a11, a12, a21, a22);
}

/* R-bin index at alpha + epsilon (nudged to land strictly inside the bin). */
static inline int ir_from_alpha(double const ds[3], double r_min, double dr, double px, double py, double alpha) {
    double ae0; /* alpha nudged by eps/|dx| in x */
    double ae1; /* alpha nudged by eps/|dy| in y */
    double x;   /* x coordinate at ae0 */
    double y;   /* y coordinate at ae1 */
    double r;   /* radial distance at nudged alpha */

    ae0 = alpha + CYL_EPS / fabs(ds[0]);
    ae1 = alpha + CYL_EPS / fabs(ds[1]);
    x = px + ae0 * ds[0];
    y = py + ae1 * ds[1];
    r = sqrt(x * x + y * y);
    return (int) floor((r - r_min) / dr);
}

/* Z-bin index at alpha + epsilon. */
static inline int iz_from_alpha(double const ds[3], double pz, double z_min, double dz, double alpha) {
    double ae; /* alpha nudged by eps/|dz| */
    double z;  /* z coordinate at ae */

    ae = alpha + CYL_EPS / fabs(ds[2]);
    z = pz + ae * ds[2];
    return (int) floor((z - z_min) / dz);
}

/* ---- main traversal ------------------------------------------------------ */

int osh_raytrace_cyl_traverse(struct osh_raytrace_grid const *grid,
                              double const p[3],
                              double const v[3],
                              double ds,
                              struct osh_voxel_crossing *crossings,
                              size_t *n_out) {
    double r_min, r_max;   /* radial grid bounds [cm] */
    double z_min, z_max;   /* axial  grid bounds [cm] */
    double dr, dz;         /* bin widths [cm] */
    size_t nr, nz;         /* bin counts */
    double ds_vec[3];      /* full displacement vector: v[i]*ds */
    double ga, gb, gc;     /* coefficients of r²(α) = ga·α²+gb·α+gc */
    double r2_p, r2_q;     /* r² at step start and end */
    double g2_min, g2_max; /* r_min² and r_max² */
    double armin, armax;   /* alpha entry/exit for the outer cylinder */
    double azmin, azmax;   /* alpha entry/exit for the Z slab */
    double az1, az2;       /* raw Z-boundary alphas before min/max */
    double amin, amax;     /* effective alpha interval to traverse */
    double a1, a2;         /* scratch roots from slv2 */
    size_t ir, iz;         /* current R-bin and Z-bin indices */
    double aur;            /* alpha of next R-shell boundary */
    double auz;            /* alpha of next Z-plane boundary */
    double ac;             /* alpha at start of current segment */
    double seg;            /* path length through current voxel [cm] */
    size_t n;              /* number of crossings accumulated */
    int i;                 /* loop counter */
    int in_hole;           /* 1 while r < r_min (inside hollow hole) */
    int ir_int, iz_int;    /* signed bin indices used before size_t cast */
    int ir_new, iz_new;    /* updated bin indices after each crossing */

    *n_out = 0;

    r_min = grid->origin[0];
    dr = grid->spacing[0];
    nr = grid->n[0];
    z_min = grid->origin[2];
    dz = grid->spacing[2];
    nz = grid->n[2];
    r_max = r_min + dr * (double) nr;
    z_max = z_min + dz * (double) nz;

    /* Full displacement vector (alpha=1 corresponds to full step ds). */
    for (i = 0; i < 3; ++i)
        ds_vec[i] = v[i] * ds;

    /* --- Quick Z-bounds rejection --- */
    if (p[2] < z_min && p[2] + ds_vec[2] < z_min)
        return 0;
    if (p[2] > z_max && p[2] + ds_vec[2] > z_max)
        return 0;

    /* --- Quick R-bounds rejection (in R²) --- */
    r2_p = p[0] * p[0] + p[1] * p[1];
    r2_q = (p[0] + ds_vec[0]) * (p[0] + ds_vec[0]) + (p[1] + ds_vec[1]) * (p[1] + ds_vec[1]);
    g2_min = r_min * r_min;
    g2_max = r_max * r_max;

    /* Both endpoints inside the hollow hole → miss. */
    if (r2_p < g2_min && r2_q < g2_min)
        return 0;

    /* Parameters for quadratic circle intersection: r²(α) = ga·α²+gb·α+gc */
    ga = ds_vec[0] * ds_vec[0] + ds_vec[1] * ds_vec[1];
    gb = 2.0 * (p[0] * ds_vec[0] + p[1] * ds_vec[1]);
    gc = p[0] * p[0] + p[1] * p[1];

    /* Both endpoints outside r_max: check if the step clips r_max at all. */
    if (r2_p > g2_max && r2_q > g2_max) {
        if (det2(ga, gb, gc, g2_max) <= 0.0)
            return 0;
        slv2(ga, gb, gc, g2_max, &a1, &a2);
        if (CYL_MIN(a1, a2) > 1.0)
            return 0;
    }

    /* Avoid divide-by-zero for purely-radial or purely-axial steps. */
    for (i = 0; i < 3; ++i) {
        if (fabs(ds_vec[i]) < CYL_EPS)
            ds_vec[i] = (ds_vec[i] >= 0.0) ? CYL_EPS : -CYL_EPS;
    }

    /* Recompute after possible clamping. */
    ga = ds_vec[0] * ds_vec[0] + ds_vec[1] * ds_vec[1];
    gb = 2.0 * (p[0] * ds_vec[0] + p[1] * ds_vec[1]);

    /* --- Alpha intervals for r_max (outer) and Z bounds --- */
    slv2(ga, gb, gc, g2_max, &a1, &a2);
    armin = CYL_MIN(a1, a2);
    armax = CYL_MAX(a1, a2);

    if (ga < CYL_EPS * CYL_EPS) { /* step parallel to Z axis */
        armin = 0.0;
        armax = 1.0;
    }
    if (CYL_MIN(armin, armax) > 1.0)
        return 0;
    if (CYL_MAX(armin, armax) < 0.0)
        return 0;

    az1 = (z_min - p[2]) / ds_vec[2];
    az2 = (z_max - p[2]) / ds_vec[2];
    azmin = CYL_MIN(az1, az2);
    azmax = CYL_MAX(az1, az2);

    amin = CYL_MAX3(0.0, armin, azmin);
    amax = CYL_MIN3(1.0, armax, azmax);
    if (amin >= amax)
        return 0;

    /* --- Starting bin indices --- */
    ir_int = ir_from_alpha(ds_vec, r_min, dr, p[0], p[1], amin);
    iz_int = iz_from_alpha(ds_vec, p[2], z_min, dz, amin);
    in_hole = (ir_int < 0) ? 1 : 0;
    if (ir_int < 0)
        ir_int = 0;
    if (iz_int < 0)
        iz_int = 0;
    if ((size_t) ir_int >= nr)
        ir_int = (int) nr - 1;
    if ((size_t) iz_int >= nz)
        iz_int = (int) nz - 1;
    ir = (size_t) ir_int;
    iz = (size_t) iz_int;

    auz = next_alpha_z(ds_vec, p[2], z_min, dz, iz);
    aur = next_alpha_r(ga, gb, gc, r_min, dr, ir, amin);

    /* --- Main traversal loop --- */
    n = 0;
    ac = amin;

    while (CYL_MIN(aur, auz) < amax) {
        /* Overflow guard: crossing buffer is pre-sized to 2*nr+nz by the caller. */
        if (n >= 2u * nr + nz)
            break;

        if (auz < aur) {
            /* Crossing a Z plane. */
            seg = (auz - ac) * ds;
            if (seg > 0.0 && !in_hole && valid_index(ir, iz, nr, nz)) {
                crossings[n].idx = flatten_idx(ir, iz, nr);
                crossings[n].path_len = seg;
                crossings[n].vol_inv = 0.0; /* filled by scoring layer */
                ++n;
            }
            ac = auz;
            iz_new = iz_from_alpha(ds_vec, p[2], z_min, dz, ac);
            if (iz_new < 0)
                iz_new = 0;
            iz = (iz_new < (int) nz) ? (size_t) iz_new : nz - 1u;
            auz = next_alpha_z(ds_vec, p[2], z_min, dz, iz);
        } else {
            /* Crossing an R shell (also handles the aur==auz tie). */
            seg = (aur - ac) * ds;
            if (seg > 0.0 && !in_hole && valid_index(ir, iz, nr, nz)) {
                crossings[n].idx = flatten_idx(ir, iz, nr);
                crossings[n].path_len = seg;
                crossings[n].vol_inv = 0.0;
                ++n;
            }
            ac = aur;
            ir_new = ir_from_alpha(ds_vec, r_min, dr, p[0], p[1], ac);
            in_hole = (ir_new < 0) ? 1 : 0;
            if (ir_new < 0)
                ir_new = 0;
            ir = (ir_new < (int) nr) ? (size_t) ir_new : nr - 1u;
            aur = next_alpha_r(ga, gb, gc, r_min, dr, ir, ac);
        }
    }

    /* Tail segment from last crossing to amax. */
    seg = (amax - ac) * ds;
    if (n < 2u * nr + nz && seg > 0.0 && !in_hole && valid_index(ir, iz, nr, nz)) {
        crossings[n].idx = flatten_idx(ir, iz, nr);
        crossings[n].path_len = seg;
        crossings[n].vol_inv = 0.0;
        ++n;
    }

    *n_out = n;
    return (n > 0u) ? 1 : 0;
}
