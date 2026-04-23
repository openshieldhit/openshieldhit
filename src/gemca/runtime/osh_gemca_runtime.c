#include "gemca/runtime/osh_gemca_runtime.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_coord.h"
#include "common/osh_diag.h"
#include "common/osh_ray.h"
#include "common/osh_vect.h"
#include "gemca/osh_gemca2.h"
#include "gemca/osh_gemca2_defines.h"
#include "gemca/runtime/osh_gemca_runtime_voxel.h"

/* ---- Local constants ----------------------------------------------------- */

/*
 * Convenience macro for selecting the smaller positive value of two doubles.
 * Identical in purpose to the MIN macro in osh_gemca2_dist.c but scoped to
 * this translation unit.
 */
#define _RT_MIN(a, b) (((a) < (b)) ? (a) : (b))
#define OSH_GEMCA_RT_BODY_BATCH_CHUNK 256

/* ---- Local types --------------------------------------------------------- */

/*
 * Stack frame used by the RPN distance evaluator.
 * Each PUSH instruction pushes one frame; binary operators pop two and push one.
 */
struct dist_frame {
    double dist;
    int is_inside;
};

/* ---- Static function prototypes ------------------------------------------ */

/* Setup helpers */
static enum osh_status setup_surfaces(struct osh_gemca_prepared const *wg, struct osh_gemca_runtime *rt);
static enum osh_status setup_bodies(struct osh_gemca_prepared const *wg, struct osh_gemca_runtime *rt);
static enum osh_status
setup_zones(struct osh_gemca_prepared const *wg, struct osh_diag_sink const *diag, struct osh_gemca_runtime *rt);
static enum osh_status setup_hu_luts(struct osh_gemca_prepared const *wg, struct osh_gemca_runtime *rt);
static int detect_zone_batch_dispatch(void);

/* Zone compilation */
static enum osh_status compile_zone(struct zone const *z,
                                    struct osh_gemca_prepared const *wg,
                                    struct osh_diag_sink const *diag,
                                    struct gemca_rt_zone *zrt);
static void compile_node(struct cgnode const *node,
                         struct osh_gemca_prepared const *wg,
                         struct osh_diag_sink const *diag,
                         struct gemca_rt_insn *insns,
                         int *ninsns);
static int count_leaves(struct cgnode const *node);
static struct body const *find_guard_body(struct cgnode const *node);
static int find_body_index(struct osh_gemca_prepared const *wg, struct body const *b);

/* RPN evaluators — inline so the compiler can fold body evaluations across the
 * guard/push boundary and into the zone-iteration outer loop. */
static inline int
eval_membership(struct osh_gemca_runtime const *rt, struct gemca_rt_zone const *z, struct ray const *r);
static void eval_membership_batch_active(struct osh_gemca_runtime const *rt,
                                         struct gemca_rt_zone const *z,
                                         double const *x,
                                         double const *y,
                                         double const *zpos,
                                         double const *ux,
                                         double const *uy,
                                         double const *uz,
                                         size_t const *candidate_idx,
                                         size_t n_candidates,
                                         int *inside_out);
static inline double
eval_distance(struct osh_gemca_runtime const *rt, struct gemca_rt_zone const *z, struct ray const *r, int *is_inside);

/* Body evaluators — inline for the same reason. */
static inline int in_body_rt(struct osh_gemca_runtime const *rt, int body_idx, struct ray const *r);
static inline double dist_body_rt(struct osh_gemca_runtime const *rt, int body_idx, struct ray const *r);

/* Ray transform — inline: called once per body check; the switch is cheap when inlined. */
static inline enum osh_status transform_to_local_rt(struct gemca_rt_body const *b, struct ray const *r, struct ray *tr);
static inline enum osh_status transform_to_local_batch_rt(struct gemca_rt_body const *b,
                                                          double const *x,
                                                          double const *y,
                                                          double const *z,
                                                          double const *ux,
                                                          double const *uy,
                                                          double const *uz,
                                                          size_t n,
                                                          double *tx,
                                                          double *ty,
                                                          double *tz,
                                                          double *tux,
                                                          double *tuy,
                                                          double *tuz);
static void check_surface_batch_indexed_rt(struct gemca_rt_surface const *sf,
                                           double const *x,
                                           double const *y,
                                           double const *z,
                                           double const *ux,
                                           double const *uy,
                                           double const *uz,
                                           size_t const *indices,
                                           size_t n,
                                           int *inside_out);
static void check_body_batch_indexed_rt(struct osh_gemca_runtime const *rt,
                                        size_t body_idx,
                                        double const *x,
                                        double const *y,
                                        double const *z,
                                        double const *ux,
                                        double const *uy,
                                        double const *uz,
                                        size_t const *indices,
                                        size_t n,
                                        int *inside_out);

/* Surface evaluators */
static inline int _check_surface_components_rt(
    struct gemca_rt_surface const *sf, double px, double py, double pz, double ux, double uy, double uz);
static inline int _check_surface_rt(struct gemca_rt_surface const *sf, struct ray const *r);
static inline double _dist_surface_rt(struct gemca_rt_surface const *sf, struct ray const *r);

/* Math helpers (pure, no state) */
static inline double _dist_sphere_rt(double r2, struct ray const *r);
static inline double _dist_cyl_rt(double r2, struct ray const *r);
static inline double _dist_elipcyl_rt(double ra2, double rb2, struct ray const *r);
static inline double _dist_cone_rt(double ra2, double rb2, struct ray const *r);
static inline double _dist_ellipsoid_rt(double ra2, double rb2, double rc2, struct ray const *r);
static inline double _dist_plane_xyz_rt(int axis, struct gemca_rt_surface const *sf, struct ray const *r);
static inline double _dist_plane_rt(struct gemca_rt_surface const *sf, struct ray const *r);
static inline double _quad_solver(double a, double b, double c);
static inline double _minpos(double a, double b);

/* ---- Public API: Lifecycle ----------------------------------------------- */

/**
 * @brief Compile a cold gemca workspace into a flat runtime representation.
 *
 * @details
 * Delegates to three setup helpers in order:
 *   1. setup_surfaces  — flatten all body surfaces into a contiguous array.
 *   2. setup_bodies    — flatten body metadata and record surface offsets.
 *   3. setup_zones     — compile each zone's cgnode AST into an RPN array.
 *
 * The workspace pointer is stored as a non-owning reference for diagnostic use.
 *
 * @param[in]  wg  Fully loaded and material-resolved cold workspace.
 * @param[out] rt  Zero-initialised runtime struct to populate.
 *
 * @returns OSH_OK on success, OSH_E* on failure.
 */
enum osh_status
osh_gemca_compile(struct osh_gemca_prepared const *wg, struct osh_diag_sink const *diag, struct osh_gemca_runtime *rt) {
    enum osh_status rc;

    if (!wg || !rt) {
        return OSH_EINVAL;
    }

    rt->workspace = wg;
    rt->zone_batch_dispatch = detect_zone_batch_dispatch();

    rc = setup_surfaces(wg, rt);
    if (rc != OSH_OK) {
        osh_gemca_runtime_free(rt);
        return rc;
    }

    rc = setup_bodies(wg, rt);
    if (rc != OSH_OK) {
        osh_gemca_runtime_free(rt);
        return rc;
    }

    rc = setup_zones(wg, diag, rt);
    if (rc != OSH_OK) {
        osh_gemca_runtime_free(rt);
        return rc;
    }

    rc = setup_hu_luts(wg, rt);
    if (rc != OSH_OK) {
        osh_gemca_runtime_free(rt);
        return rc;
    }

    return OSH_OK;
}

char const *osh_gemca_runtime_zone_batch_dispatch_name(struct osh_gemca_runtime const *rt) {
    if (!rt) {
        return "scalar";
    }

    switch (rt->zone_batch_dispatch) {
    case OSH_GEMCA_ZONE_BATCH_DISPATCH_AVX2:
        return "AVX2+FMA (4-wide SIMD)";
    case OSH_GEMCA_ZONE_BATCH_DISPATCH_SCALAR_NOCPU:
        return "scalar (AVX2 built-in but CPU lacks support)";
    case OSH_GEMCA_ZONE_BATCH_DISPATCH_SCALAR:
    default:
        return "scalar";
    }
}

/**
 * @brief Release all allocations owned by a gemca runtime.
 *
 * @details
 * Frees per-zone insns[] buffers first, then zones[], bodies[], surfaces[].
 * The cold workspace pointer is not freed.  Safe to call on a zero-initialised
 * or partially initialised struct.
 *
 * @param[in] rt  Runtime to release; may be NULL.
 */
void osh_gemca_runtime_free(struct osh_gemca_runtime *rt) {
    size_t i;

    if (!rt) {
        return;
    }

    if (rt->zones) {
        for (i = 0; i < rt->nzones; i++) {
            free(rt->zones[i].insns);
            rt->zones[i].insns = NULL;
        }
        free(rt->zones);
        rt->zones = NULL;
    }

    free(rt->bodies);
    rt->bodies = NULL;

    free(rt->surfaces);
    rt->surfaces = NULL;

    free(rt->hu_bin_lut);
    rt->hu_bin_lut = NULL;

    free(rt->hu_rho_lut);
    rt->hu_rho_lut = NULL;

    rt->nsurfaces = 0;
    rt->nbodies = 0;
    rt->nzones = 0;
}

/* ---- Public API: Hot query ---------------------------------------------- */

/**
 * @brief Return the dense zone index for the zone containing a ray.
 *
 * @details
 * Iterates zones in order; evaluates each zone's RPN expression via
 * eval_membership().  GUARD_BODY instructions at the head of the array provide
 * O(1) rejection for most zones.
 *
 * @param[in] rt  Compiled gemca runtime.
 * @param[in] r   Ray to locate.
 *
 * @returns 0-based zone index, or OSH_GEMCA_ZONE_INDEX_INVALID.
 */
size_t osh_gemca_runtime_get_zone(struct osh_gemca_runtime const *rt, struct ray const *r) {
    size_t i;

    if (!rt || !r) {
        return OSH_GEMCA_ZONE_INDEX_INVALID;
    }

    for (i = 0; i < rt->nzones; i++) {
        if (eval_membership(rt, &rt->zones[i], r)) {
            return i;
        }
    }

    return OSH_GEMCA_ZONE_INDEX_INVALID;
}

/**
 * @brief Return the total distance a ray travels inside zone `zone_idx`.
 *
 * @details
 * Normalises the ray direction on a local copy, then advances the copy through
 * the zone one boundary at a time, accumulating path length, until the ray
 * leaves the zone (is_inside == 0).
 *
 * @param[in] rt        Compiled gemca runtime.
 * @param[in] zone_idx  Zone index from osh_gemca_runtime_get_zone().
 * @param[in] r         Ray (position and direction; caller's struct is unmodified).
 *
 * @returns Total path length inside the zone.
 */
double osh_gemca_runtime_get_distance(struct osh_gemca_runtime const *rt, size_t zone_idx, struct ray const *r) {
    struct ray rr;
    double d;
    double total;
    int is_inside;
    int i;

    if (!rt || !r || zone_idx >= rt->nzones) {
        return 0.0;
    }

    /* Work on a local copy so the caller's ray is unchanged. */
    for (i = 0; i < 3; i++) {
        rr.p[i] = r->p[i];
        rr.cp[i] = r->cp[i];
    }
    rr.system = r->system;
    osh_vect_norm(rr.cp);

    total = 0.0;

    while (1) {
        d = eval_distance(rt, &rt->zones[zone_idx], &rr, &is_inside);
        if (!is_inside) {
            break;
        }
        if (d < OSH_GEMCA_STEPLIM) {
            d = OSH_GEMCA_STEPLIM;
        }
        total += d;
        /* Advance ray position along normalised direction. */
        for (i = 0; i < 3; i++) {
            rr.p[i] += rr.cp[i] * d;
        }
    }

    return total;
}

/* ---- Batch query API ----------------------------------------------------- */

static void check_surface_batch_indexed_rt(struct gemca_rt_surface const *sf,
                                           double const *x,
                                           double const *y,
                                           double const *z,
                                           double const *ux,
                                           double const *uy,
                                           double const *uz,
                                           size_t const *indices,
                                           size_t n,
                                           int *inside_out) {
    double d;
    double dot;
    size_t lane;
    size_t i;
    double p0;
    double p1;
    double p2;
    double p3;

    if (!sf || !x || !y || !z || !ux || !uy || !uz || !indices || !inside_out) {
        return;
    }

    p0 = sf->p[0];
    p1 = sf->p[1];
    p2 = sf->p[2];
    p3 = sf->p[3];

    switch (sf->type) {
    case OSH_GEMCA_SURF_SPHERE:
        for (lane = 0; lane < n; ++lane) {
            i = indices[lane];
            d = (x[i] * x[i]) + (y[i] * y[i]) + (z[i] * z[i]) - p0;
            if (d > OSH_GEMCA_SMALL) {
                inside_out[lane] = 0;
            } else if (d < -OSH_GEMCA_SMALL) {
                inside_out[lane] = 1;
            } else {
                inside_out[lane] = (((x[i] * ux[i]) + (y[i] * uy[i]) + (z[i] * uz[i])) < 0.0) ? 1 : 0;
            }
        }
        break;

    case OSH_GEMCA_SURF_ELLIPSOID:
        for (lane = 0; lane < n; ++lane) {
            i = indices[lane];
            d = (x[i] * x[i]) / p0 + (y[i] * y[i]) / p1 + (z[i] * z[i]) / p2;
            if (d > 1.0 + OSH_GEMCA_SMALL) {
                inside_out[lane] = 0;
            } else if (d < 1.0 - OSH_GEMCA_SMALL) {
                inside_out[lane] = 1;
            } else {
                dot = (x[i] / p0) * ux[i] + (y[i] / p1) * uy[i] + (z[i] / p2) * uz[i];
                inside_out[lane] = (dot < 0.0) ? 1 : 0;
            }
        }
        break;

    case OSH_GEMCA_SURF_CYLZ:
        for (lane = 0; lane < n; ++lane) {
            i = indices[lane];
            d = (x[i] * x[i]) + (y[i] * y[i]) - p0;
            if (d > OSH_GEMCA_SMALL) {
                inside_out[lane] = 0;
            } else if (d < -OSH_GEMCA_SMALL) {
                inside_out[lane] = 1;
            } else {
                dot = (x[i] * ux[i]) + (y[i] * uy[i]);
                inside_out[lane] = (dot < 0.0) ? 1 : 0;
            }
        }
        break;

    case OSH_GEMCA_SURF_ELLZ:
        for (lane = 0; lane < n; ++lane) {
            i = indices[lane];
            d = (x[i] * x[i] / p0) + (y[i] * y[i] / p1) - 1.0;
            if (d > OSH_GEMCA_SMALL) {
                inside_out[lane] = 0;
            } else if (d < -OSH_GEMCA_SMALL) {
                inside_out[lane] = 1;
            } else {
                dot = (x[i] / p0) * ux[i] + (y[i] / p1) * uy[i];
                inside_out[lane] = (dot < 0.0) ? 1 : 0;
            }
        }
        break;

    case OSH_GEMCA_SURF_CONE:
        for (lane = 0; lane < n; ++lane) {
            i = indices[lane];
            d = (x[i] * x[i]) + (y[i] * y[i]) - p1 * (z[i] * z[i]);
            if (d > OSH_GEMCA_SMALL) {
                inside_out[lane] = 0;
            } else if (d < -OSH_GEMCA_SMALL) {
                inside_out[lane] = 1;
            } else {
                dot = (x[i] * ux[i]) + (y[i] * uy[i]) - (p1 * z[i] * uz[i]);
                inside_out[lane] = (dot < 0.0) ? 1 : 0;
            }
        }
        break;

    case OSH_GEMCA_SURF_PLANEX:
        for (lane = 0; lane < n; ++lane) {
            i = indices[lane];
            d = p0 * x[i] + p1;
            if (d > OSH_GEMCA_SMALL) {
                inside_out[lane] = 0;
            } else if (d < -OSH_GEMCA_SMALL) {
                inside_out[lane] = 1;
            } else {
                inside_out[lane] = (p0 * ux[i] > 0.0) ? 0 : 1;
            }
        }
        break;

    case OSH_GEMCA_SURF_PLANEY:
        for (lane = 0; lane < n; ++lane) {
            i = indices[lane];
            d = p0 * y[i] + p1;
            if (d > OSH_GEMCA_SMALL) {
                inside_out[lane] = 0;
            } else if (d < -OSH_GEMCA_SMALL) {
                inside_out[lane] = 1;
            } else {
                inside_out[lane] = (p0 * uy[i] > 0.0) ? 0 : 1;
            }
        }
        break;

    case OSH_GEMCA_SURF_PLANEZ:
        for (lane = 0; lane < n; ++lane) {
            i = indices[lane];
            d = p0 * z[i] + p1;
            if (d > OSH_GEMCA_SMALL) {
                inside_out[lane] = 0;
            } else if (d < -OSH_GEMCA_SMALL) {
                inside_out[lane] = 1;
            } else {
                inside_out[lane] = (p0 * uz[i] > 0.0) ? 0 : 1;
            }
        }
        break;

    case OSH_GEMCA_SURF_PLANE:
        for (lane = 0; lane < n; ++lane) {
            i = indices[lane];
            d = (p0 * x[i]) + (p1 * y[i]) + (p2 * z[i]) + p3;
            if (d > OSH_GEMCA_SMALL) {
                inside_out[lane] = 0;
            } else if (d < -OSH_GEMCA_SMALL) {
                inside_out[lane] = 1;
            } else {
                inside_out[lane] = ((p0 * ux[i]) + (p1 * uy[i]) + (p2 * uz[i]) > 0.0) ? 0 : 1;
            }
        }
        break;

    default:
        for (lane = 0; lane < n; ++lane) {
            i = indices[lane];
            inside_out[lane] = _check_surface_components_rt(sf, x[i], y[i], z[i], ux[i], uy[i], uz[i]);
        }
        break;
    }
}

static void check_body_batch_indexed_rt(struct osh_gemca_runtime const *rt,
                                        size_t body_idx,
                                        double const *x,
                                        double const *y,
                                        double const *z,
                                        double const *ux,
                                        double const *uy,
                                        double const *uz,
                                        size_t const *indices,
                                        size_t n,
                                        int *inside_out) {
    struct gemca_rt_body const *b;
    size_t lane;
    size_t is;
    size_t local_idx[OSH_GEMCA_RT_BODY_BATCH_CHUNK];
    double tx[OSH_GEMCA_RT_BODY_BATCH_CHUNK];
    double ty[OSH_GEMCA_RT_BODY_BATCH_CHUNK];
    double tz[OSH_GEMCA_RT_BODY_BATCH_CHUNK];
    double tux[OSH_GEMCA_RT_BODY_BATCH_CHUNK];
    double tuy[OSH_GEMCA_RT_BODY_BATCH_CHUNK];
    double tuz[OSH_GEMCA_RT_BODY_BATCH_CHUNK];
    int surface_inside[OSH_GEMCA_RT_BODY_BATCH_CHUNK];
    double const *lx;
    double const *ly;
    double const *lz;
    double const *lux;
    double const *luy;
    double const *luz;

    if (!x || !y || !z || !ux || !uy || !uz || !indices || !inside_out) {
        return;
    }

    if (!rt || body_idx >= rt->nbodies || n > (size_t) OSH_GEMCA_RT_BODY_BATCH_CHUNK) {
        for (lane = 0; lane < n; ++lane) {
            inside_out[lane] = 0;
        }
        return;
    }

    b = &rt->bodies[body_idx];
    lx = x;
    ly = y;
    lz = z;
    lux = ux;
    luy = uy;
    luz = uz;

    switch (b->coord) {
    case OSH_COORD_UNIVERSE:
        for (lane = 0; lane < n; ++lane) {
            inside_out[lane] = 1;
            local_idx[lane] = indices[lane];
        }
        break;

    case OSH_COORD_BCALIGN:
        for (lane = 0; lane < n; ++lane) {
            size_t i = indices[lane];
            tx[lane] = x[i] + b->t[3];
            ty[lane] = y[i] + b->t[7];
            tz[lane] = z[i] + b->t[11];
            tux[lane] = ux[i];
            tuy[lane] = uy[i];
            tuz[lane] = uz[i];
            inside_out[lane] = 1;
            local_idx[lane] = lane;
        }
        lx = tx;
        ly = ty;
        lz = tz;
        lux = tux;
        luy = tuy;
        luz = tuz;
        break;

    case OSH_COORD_BZALIGN:
        for (lane = 0; lane < n; ++lane) {
            size_t i = indices[lane];
            tx[lane] = x[i] * b->t[0] + y[i] * b->t[1] + z[i] * b->t[2] - b->t[3];
            ty[lane] = x[i] * b->t[4] + y[i] * b->t[5] + z[i] * b->t[6] - b->t[7];
            tz[lane] = x[i] * b->t[8] + y[i] * b->t[9] + z[i] * b->t[10] - b->t[11];
            tux[lane] = ux[i] * b->t[0] + uy[i] * b->t[1] + uz[i] * b->t[2];
            tuy[lane] = ux[i] * b->t[4] + uy[i] * b->t[5] + uz[i] * b->t[6];
            tuz[lane] = ux[i] * b->t[8] + uy[i] * b->t[9] + uz[i] * b->t[10];
            inside_out[lane] = 1;
            local_idx[lane] = lane;
        }
        lx = tx;
        ly = ty;
        lz = tz;
        lux = tux;
        luy = tuy;
        luz = tuz;
        break;

    default:
        for (lane = 0; lane < n; ++lane) {
            inside_out[lane] = 0;
        }
        return;
    }

    for (is = 0; is < (size_t) b->nsurfs; ++is) {
        int any_inside;

        check_surface_batch_indexed_rt(
            &rt->surfaces[b->surf_begin + is], lx, ly, lz, lux, luy, luz, local_idx, n, surface_inside);

        any_inside = 0;
        for (lane = 0; lane < n; ++lane) {
            inside_out[lane] = inside_out[lane] && surface_inside[lane];
            any_inside |= inside_out[lane];
        }

        if (!any_inside) {
            break;
        }
    }
}

/**
 * @brief Surface-membership query for a batch of @p n rays (SoA inputs).
 *
 * @details
 * Thin batched wrapper over the scalar per-primitive formulas shared with
 * _check_surface_rt().  Inputs must already be expressed in the surface's
 * local body coordinates.
 */
void osh_gemca_runtime_check_surface_batch(struct gemca_rt_surface const *sf,
                                           double const *x,
                                           double const *y,
                                           double const *z,
                                           double const *ux,
                                           double const *uy,
                                           double const *uz,
                                           size_t n,
                                           int *inside_out) {
    double d;
    double dot;
    size_t i;
    double p0;
    double p1;
    double p2;
    double p3;

    if (!sf || !x || !y || !z || !ux || !uy || !uz || !inside_out) {
        return;
    }

    p0 = sf->p[0];
    p1 = sf->p[1];
    p2 = sf->p[2];
    p3 = sf->p[3];

    switch (sf->type) {
    case OSH_GEMCA_SURF_SPHERE:
        for (i = 0; i < n; ++i) {
            d = (x[i] * x[i]) + (y[i] * y[i]) + (z[i] * z[i]) - p0;
            if (d > OSH_GEMCA_SMALL) {
                inside_out[i] = 0;
            } else if (d < -OSH_GEMCA_SMALL) {
                inside_out[i] = 1;
            } else {
                inside_out[i] = (((x[i] * ux[i]) + (y[i] * uy[i]) + (z[i] * uz[i])) < 0.0) ? 1 : 0;
            }
        }
        break;

    case OSH_GEMCA_SURF_ELLIPSOID:
        for (i = 0; i < n; ++i) {
            d = (x[i] * x[i]) / p0 + (y[i] * y[i]) / p1 + (z[i] * z[i]) / p2;
            if (d > 1.0 + OSH_GEMCA_SMALL) {
                inside_out[i] = 0;
            } else if (d < 1.0 - OSH_GEMCA_SMALL) {
                inside_out[i] = 1;
            } else {
                dot = (x[i] / p0) * ux[i] + (y[i] / p1) * uy[i] + (z[i] / p2) * uz[i];
                inside_out[i] = (dot < 0.0) ? 1 : 0;
            }
        }
        break;

    case OSH_GEMCA_SURF_CYLZ:
        for (i = 0; i < n; ++i) {
            d = (x[i] * x[i]) + (y[i] * y[i]) - p0;
            if (d > OSH_GEMCA_SMALL) {
                inside_out[i] = 0;
            } else if (d < -OSH_GEMCA_SMALL) {
                inside_out[i] = 1;
            } else {
                dot = (x[i] * ux[i]) + (y[i] * uy[i]);
                inside_out[i] = (dot < 0.0) ? 1 : 0;
            }
        }
        break;

    case OSH_GEMCA_SURF_ELLZ:
        for (i = 0; i < n; ++i) {
            d = (x[i] * x[i] / p0) + (y[i] * y[i] / p1) - 1.0;
            if (d > OSH_GEMCA_SMALL) {
                inside_out[i] = 0;
            } else if (d < -OSH_GEMCA_SMALL) {
                inside_out[i] = 1;
            } else {
                dot = (x[i] / p0) * ux[i] + (y[i] / p1) * uy[i];
                inside_out[i] = (dot < 0.0) ? 1 : 0;
            }
        }
        break;

    case OSH_GEMCA_SURF_CONE:
        for (i = 0; i < n; ++i) {
            d = (x[i] * x[i]) + (y[i] * y[i]) - p1 * (z[i] * z[i]);
            if (d > OSH_GEMCA_SMALL) {
                inside_out[i] = 0;
            } else if (d < -OSH_GEMCA_SMALL) {
                inside_out[i] = 1;
            } else {
                dot = (x[i] * ux[i]) + (y[i] * uy[i]) - (p1 * z[i] * uz[i]);
                inside_out[i] = (dot < 0.0) ? 1 : 0;
            }
        }
        break;

    case OSH_GEMCA_SURF_PLANEX:
        for (i = 0; i < n; ++i) {
            d = p0 * x[i] + p1;
            if (d > OSH_GEMCA_SMALL) {
                inside_out[i] = 0;
            } else if (d < -OSH_GEMCA_SMALL) {
                inside_out[i] = 1;
            } else {
                inside_out[i] = (p0 * ux[i] > 0.0) ? 0 : 1;
            }
        }
        break;

    case OSH_GEMCA_SURF_PLANEY:
        for (i = 0; i < n; ++i) {
            d = p0 * y[i] + p1;
            if (d > OSH_GEMCA_SMALL) {
                inside_out[i] = 0;
            } else if (d < -OSH_GEMCA_SMALL) {
                inside_out[i] = 1;
            } else {
                inside_out[i] = (p0 * uy[i] > 0.0) ? 0 : 1;
            }
        }
        break;

    case OSH_GEMCA_SURF_PLANEZ:
        for (i = 0; i < n; ++i) {
            d = p0 * z[i] + p1;
            if (d > OSH_GEMCA_SMALL) {
                inside_out[i] = 0;
            } else if (d < -OSH_GEMCA_SMALL) {
                inside_out[i] = 1;
            } else {
                inside_out[i] = (p0 * uz[i] > 0.0) ? 0 : 1;
            }
        }
        break;

    case OSH_GEMCA_SURF_PLANE:
        for (i = 0; i < n; ++i) {
            d = (p0 * x[i]) + (p1 * y[i]) + (p2 * z[i]) + p3;
            if (d > OSH_GEMCA_SMALL) {
                inside_out[i] = 0;
            } else if (d < -OSH_GEMCA_SMALL) {
                inside_out[i] = 1;
            } else {
                inside_out[i] = ((p0 * ux[i]) + (p1 * uy[i]) + (p2 * uz[i]) > 0.0) ? 0 : 1;
            }
        }
        break;

    default:
        for (i = 0; i < n; ++i) {
            inside_out[i] = _check_surface_components_rt(sf, x[i], y[i], z[i], ux[i], uy[i], uz[i]);
        }
        break;
    }
}

/**
 * @brief Body-membership query for a batch of @p n rays (SoA inputs).
 *
 * @details
 * Processes the batch in fixed-size chunks to keep the temporary transformed
 * coordinates on the stack.  Each chunk is transformed once into the body's
 * local coordinate system, then every surface predicate is applied across the
 * chunk via osh_gemca_runtime_check_surface_batch().
 */
void osh_gemca_runtime_check_body_batch(struct osh_gemca_runtime const *rt,
                                        size_t body_idx,
                                        double const *x,
                                        double const *y,
                                        double const *z,
                                        double const *ux,
                                        double const *uy,
                                        double const *uz,
                                        size_t n,
                                        int *inside_out) {
    struct gemca_rt_body const *b;
    double const *lx;
    double const *ly;
    double const *lz;
    double const *lux;
    double const *luy;
    double const *luz;
    size_t base;
    size_t i;
    size_t is;
    size_t chunk;
    int any_inside;
    double tx[OSH_GEMCA_RT_BODY_BATCH_CHUNK];
    double ty[OSH_GEMCA_RT_BODY_BATCH_CHUNK];
    double tz[OSH_GEMCA_RT_BODY_BATCH_CHUNK];
    double tux[OSH_GEMCA_RT_BODY_BATCH_CHUNK];
    double tuy[OSH_GEMCA_RT_BODY_BATCH_CHUNK];
    double tuz[OSH_GEMCA_RT_BODY_BATCH_CHUNK];
    int surface_inside[OSH_GEMCA_RT_BODY_BATCH_CHUNK];

    if (!x || !y || !z || !ux || !uy || !uz || !inside_out) {
        return;
    }

    if (!rt || body_idx >= rt->nbodies) {
        for (i = 0; i < n; ++i) {
            inside_out[i] = 0;
        }
        return;
    }

    b = &rt->bodies[body_idx];

    for (base = 0; base < n; base += chunk) {
        int *inside_chunk;

        chunk = n - base;
        if (chunk > (size_t) OSH_GEMCA_RT_BODY_BATCH_CHUNK) {
            chunk = (size_t) OSH_GEMCA_RT_BODY_BATCH_CHUNK;
        }

        inside_chunk = inside_out + base;
        lx = x + base;
        ly = y + base;
        lz = z + base;
        lux = ux + base;
        luy = uy + base;
        luz = uz + base;

        switch (b->coord) {
        case OSH_COORD_UNIVERSE:
            break;

        case OSH_COORD_BCALIGN:
            for (i = 0; i < chunk; ++i) {
                tx[i] = lx[i] + b->t[3];
                ty[i] = ly[i] + b->t[7];
                tz[i] = lz[i] + b->t[11];
            }
            lx = tx;
            ly = ty;
            lz = tz;
            break;

        case OSH_COORD_BZALIGN:
            if (transform_to_local_batch_rt(b, lx, ly, lz, lux, luy, luz, chunk, tx, ty, tz, tux, tuy, tuz) != OSH_OK) {
                for (i = 0; i < chunk; ++i) {
                    inside_chunk[i] = 0;
                }
                continue;
            }
            lx = tx;
            ly = ty;
            lz = tz;
            lux = tux;
            luy = tuy;
            luz = tuz;
            break;

        default:
            for (i = 0; i < chunk; ++i) {
                inside_chunk[i] = 0;
            }
            continue;
        }

        for (i = 0; i < chunk; ++i) {
            inside_chunk[i] = 1;
        }

        for (is = 0; is < (size_t) b->nsurfs; ++is) {
            osh_gemca_runtime_check_surface_batch(
                &rt->surfaces[b->surf_begin + is], lx, ly, lz, lux, luy, luz, chunk, surface_inside);

            any_inside = 0;
            for (i = 0; i < chunk; ++i) {
                inside_chunk[i] = inside_chunk[i] && surface_inside[i];
                any_inside |= inside_chunk[i];
            }

            if (!any_inside) {
                break;
            }
        }
    }
}

/* Forward declaration for AVX2 implementation (compiled separately). */
#ifdef OSH_GEMCA_RUNTIME_HAVE_AVX2
extern void osh_gemca_runtime_get_zone_batch_avx2(struct osh_gemca_runtime const *rt,
                                                  double const *x,
                                                  double const *y,
                                                  double const *z,
                                                  double const *ux,
                                                  double const *uy,
                                                  double const *uz,
                                                  size_t n,
                                                  size_t *zone_out);
#endif

/**
 * @brief Zone lookup for a batch of @p n particles (SoA inputs).
 *
 * @details
 * Dispatches to the AVX2+FMA implementation when the CPU supports it (runtime
 * check via __builtin_cpu_supports), otherwise falls back to the scalar
 * chunked-active-list implementation.
 */
void osh_gemca_runtime_get_zone_batch(struct osh_gemca_runtime const *rt,
                                      double const *x,
                                      double const *y,
                                      double const *z,
                                      double const *ux,
                                      double const *uy,
                                      double const *uz,
                                      size_t n,
                                      size_t *zone_out) {
    size_t base;
    size_t i;
    size_t j;
    size_t chunk;
    size_t unresolved_n;
    size_t unresolved_idx[OSH_GEMCA_RT_BODY_BATCH_CHUNK];
    int inside_chunk[OSH_GEMCA_RT_BODY_BATCH_CHUNK];

    if (!rt || !x || !y || !z || !ux || !uy || !uz || !zone_out) {
        return;
    }

#ifdef OSH_GEMCA_RUNTIME_HAVE_AVX2
    if (__builtin_cpu_supports("avx2")) {
        osh_gemca_runtime_get_zone_batch_avx2(rt, x, y, z, ux, uy, uz, n, zone_out);
        return;
    }
#endif

    for (base = 0; base < n; base += chunk) {
        chunk = n - base;
        if (chunk > (size_t) OSH_GEMCA_RT_BODY_BATCH_CHUNK) {
            chunk = (size_t) OSH_GEMCA_RT_BODY_BATCH_CHUNK;
        }

        for (i = 0; i < chunk; ++i) {
            zone_out[base + i] = OSH_GEMCA_ZONE_INDEX_INVALID;
            unresolved_idx[i] = i;
        }

        unresolved_n = chunk;
        for (j = 0; j < rt->nzones && unresolved_n > 0u; ++j) {
            size_t write = 0u;
            size_t k;

            eval_membership_batch_active(rt,
                                         &rt->zones[j],
                                         x + base,
                                         y + base,
                                         z + base,
                                         ux + base,
                                         uy + base,
                                         uz + base,
                                         unresolved_idx,
                                         unresolved_n,
                                         inside_chunk);

            for (k = 0; k < unresolved_n; ++k) {
                i = unresolved_idx[k];
                if (inside_chunk[i]) {
                    zone_out[base + i] = j;
                } else {
                    unresolved_idx[write++] = i;
                }
            }
            unresolved_n = write;
        }
    }
}

static int detect_zone_batch_dispatch(void) {
#ifdef OSH_GEMCA_RUNTIME_HAVE_AVX2
    if (__builtin_cpu_supports("avx2")) {
        return OSH_GEMCA_ZONE_BATCH_DISPATCH_AVX2;
    }
    return OSH_GEMCA_ZONE_BATCH_DISPATCH_SCALAR_NOCPU;
#else
    return OSH_GEMCA_ZONE_BATCH_DISPATCH_SCALAR;
#endif
}

/**
 * @brief Boundary-distance query for a batch of @p n particles (SoA inputs).
 *
 * @details
 * Skips particles whose zone_indices entry is OSH_GEMCA_ZONE_INDEX_INVALID
 * (outside all geometry) and writes 0.0 for those slots.  All other slots
 * call eval_distance() via the same step-loop used in get_distance().
 */
void osh_gemca_runtime_get_distance_batch(struct osh_gemca_runtime const *rt,
                                          double const *x,
                                          double const *y,
                                          double const *z,
                                          double const *ux,
                                          double const *uy,
                                          double const *uz,
                                          size_t const *zone_indices,
                                          size_t n,
                                          double *dist_out) {
    struct ray rr;
    size_t i;
    int k;
    size_t zone_idx;
    double d;
    double total;
    int is_inside;

    if (!rt || !x || !y || !z || !ux || !uy || !uz || !zone_indices || !dist_out) {
        if (dist_out && n) {
            size_t _k;
            for (_k = 0; _k < n; ++_k)
                dist_out[_k] = 0.0;
        }
        return;
    }

    for (i = 0; i < n; ++i) {
        zone_idx = zone_indices[i];
        if (zone_idx == OSH_GEMCA_ZONE_INDEX_INVALID || zone_idx >= rt->nzones) {
            dist_out[i] = 0.0;
            continue;
        }

        /* Work on a local ray copy — same normalisation as scalar variant. */
        rr.p[0] = x[i];
        rr.p[1] = y[i];
        rr.p[2] = z[i];
        rr.cp[0] = ux[i];
        rr.cp[1] = uy[i];
        rr.cp[2] = uz[i];
        rr.system = OSH_COORD_UNIVERSE;
        osh_vect_norm(rr.cp);

        total = 0.0;
        while (1) {
            d = eval_distance(rt, &rt->zones[zone_idx], &rr, &is_inside);
            if (!is_inside) {
                break;
            }
            if (d < OSH_GEMCA_STEPLIM) {
                d = OSH_GEMCA_STEPLIM;
            }
            total += d;
            for (k = 0; k < 3; ++k) {
                rr.p[k] += rr.cp[k] * d;
            }
        }
        dist_out[i] = total;
    }
}

/* ---- Setup helpers ------------------------------------------------------- */

/**
 * @brief Flatten all body surfaces from the cold workspace into rt->surfaces[].
 *
 * @details
 * Counts total surfaces across all bodies, allocates a single contiguous array,
 * then copies parameters from each cold surface into a fixed-size flat entry.
 * The surface order matches the body order in wg->bodies[]: all surfaces for
 * body 0, then all for body 1, and so on.  setup_bodies() uses this ordering
 * when recording surf_begin offsets.
 *
 * @param[in]  wg  Cold workspace.
 * @param[out] rt  Runtime to populate (rt->surfaces allocated and filled).
 *
 * @returns OSH_OK or OSH_ENOMEM.
 */
static enum osh_status setup_surfaces(struct osh_gemca_prepared const *wg, struct osh_gemca_runtime *rt) {
    size_t total;
    size_t ib;
    size_t is;
    int ip;
    struct body const *b;
    struct surface const *sf;
    struct gemca_rt_surface *dst;

    total = 0;
    for (ib = 0; ib < wg->nbodies; ib++) {
        total += (size_t) wg->bodies[ib]->nsurfs;
    }

    rt->surfaces = (struct gemca_rt_surface *) calloc(total, sizeof(struct gemca_rt_surface));
    if (!rt->surfaces) {
        return OSH_ENOMEM;
    }
    rt->nsurfaces = total;

    is = 0;
    for (ib = 0; ib < wg->nbodies; ib++) {
        b = wg->bodies[ib];
        for (ip = 0; ip < b->nsurfs; ip++) {
            sf = b->surfs[ip];
            dst = &rt->surfaces[is];
            dst->type = sf->type;
            /* Copy available parameters; remaining slots stay zero-padded. */
            memcpy(dst->p, sf->p, (size_t) sf->np * sizeof(double));
            is++;
        }
    }

    return OSH_OK;
}

/**
 * @brief Compile body metadata from the cold workspace into rt->bodies[].
 *
 * @details
 * Records the transformation matrix, voxel-grid metadata, coordinate system,
 * body type, and the offset (surf_begin) into rt->surfaces[] for each body.
 * The offset is computed by accumulating surface counts in body order, which
 * matches the layout produced by setup_surfaces().
 *
 * @param[in]  wg  Cold workspace.
 * @param[out] rt  Runtime (rt->surfaces must already be populated).
 *
 * @returns OSH_OK or OSH_ENOMEM.
 */
static enum osh_status setup_hu_luts(struct osh_gemca_prepared const *wg, struct osh_gemca_runtime *rt) {
    if (!wg->hu_bin_lut || !wg->hu_rho_lut) {
        return OSH_OK;
    }

    rt->hu_bin_lut = (uint8_t *) malloc(2601u * sizeof(uint8_t));
    if (!rt->hu_bin_lut) {
        return OSH_ENOMEM;
    }
    rt->hu_rho_lut = (float *) malloc(2601u * sizeof(float));
    if (!rt->hu_rho_lut) {
        free(rt->hu_bin_lut);
        rt->hu_bin_lut = NULL;
        return OSH_ENOMEM;
    }

    memcpy(rt->hu_bin_lut, wg->hu_bin_lut, 2601u * sizeof(uint8_t));
    memcpy(rt->hu_rho_lut, wg->hu_rho_lut, 2601u * sizeof(float));

    return OSH_OK;
}

static enum osh_status setup_bodies(struct osh_gemca_prepared const *wg, struct osh_gemca_runtime *rt) {
    size_t ib;
    size_t surf_offset;
    struct body const *b;
    struct gemca_rt_body *dst;

    rt->bodies = (struct gemca_rt_body *) calloc(wg->nbodies, sizeof(struct gemca_rt_body));
    if (!rt->bodies) {
        return OSH_ENOMEM;
    }
    rt->nbodies = wg->nbodies;

    surf_offset = 0;
    for (ib = 0; ib < wg->nbodies; ib++) {
        b = wg->bodies[ib];
        dst = &rt->bodies[ib];

        memcpy(dst->t, b->t, 16 * sizeof(double));
        memcpy(&dst->ct_grid, &b->ct_grid, sizeof(dst->ct_grid));
        dst->hu = b->hu;
        dst->surf_begin = surf_offset;
        dst->nsurfs = b->nsurfs;
        dst->type = b->type;
        dst->coord = b->coord;

        surf_offset += (size_t) b->nsurfs;
    }

    return OSH_OK;
}

/**
 * @brief Compile each zone's CSG expression into a flat RPN instruction array.
 *
 * @details
 * Allocates rt->zones[] and for each zone calls compile_zone() to produce
 * the RPN instruction sequence.  Zone material indices are copied directly from
 * the cold zone struct (they must be resolved before this function is called).
 *
 * @param[in]  wg  Cold workspace.
 * @param[out] rt  Runtime (rt->surfaces and rt->bodies must be populated).
 *
 * @returns OSH_OK or OSH_ENOMEM.
 */
static enum osh_status
setup_zones(struct osh_gemca_prepared const *wg, struct osh_diag_sink const *diag, struct osh_gemca_runtime *rt) {
    size_t iz;
    enum osh_status rc;

    rt->zones = (struct gemca_rt_zone *) calloc(wg->nzones, sizeof(struct gemca_rt_zone));
    if (!rt->zones) {
        return OSH_ENOMEM;
    }
    rt->nzones = wg->nzones;

    for (iz = 0; iz < wg->nzones; iz++) {
        rc = compile_zone(wg->zones[iz], wg, diag, &rt->zones[iz]);
        if (rc != OSH_OK) {
            return rc;
        }
    }

    return OSH_OK;
}

/* ---- Zone compilation ---------------------------------------------------- */

/**
 * @brief Compile one zone's cgnode AST into a flat RPN instruction array.
 *
 * @details
 * Determines the instruction budget (2 * leaves + 1 for optional guard),
 * allocates the array, optionally prepends a GUARD_BODY instruction derived
 * from find_guard_body(), then recursively compiles the AST via compile_node().
 *
 * After compilation the array is scanned for PUSH_VOXEL_BODY to record
 * voxel_body_idx for later Jacobs dispatch.
 *
 * Guard derivation: find_guard_body() follows the leftmost branch of the tree
 * through only intersection (+) and difference (-) operators.  If that path
 * reaches a leaf, the leaf body guards the entire zone (any point outside it
 * cannot be inside the zone).  For zones rooted at a union (|), no simple
 * single-body guard exists and none is prepended.
 *
 * @param[in]  z    Cold zone to compile.
 * @param[in]  wg   Cold workspace (for body-pointer-to-index lookup).
 * @param[out] zrt  Runtime zone to populate.
 *
 * @returns OSH_OK or OSH_ENOMEM.
 */
static enum osh_status compile_zone(struct zone const *z,
                                    struct osh_gemca_prepared const *wg,
                                    struct osh_diag_sink const *diag,
                                    struct gemca_rt_zone *zrt) {
    struct body const *guard;
    int guard_idx;
    int nleaves;
    int max_insns;
    int ninsns;
    int i;

    nleaves = count_leaves(&z->node);
    /*
     * An N-leaf binary tree has N-1 interior nodes, so RPN needs exactly
     * 2*N-1 instructions.  One extra slot is reserved for an optional guard.
     */
    max_insns = 2 * nleaves + 1;

    zrt->insns = (struct gemca_rt_insn *) malloc((size_t) max_insns * sizeof(struct gemca_rt_insn));
    if (!zrt->insns) {
        return OSH_ENOMEM;
    }

    ninsns = 0;

    /*
     * Prepend a GUARD_BODY only for multi-leaf zones.
     *
     * For a single-leaf zone the only instruction is PUSH_BODY(A), which
     * already performs exactly the same test a GUARD would do.  Emitting
     * GUARD_BODY(A) PUSH_BODY(A) would evaluate body A twice per ray —
     * doubling the per-zone cost with no benefit.
     */
    if (nleaves > 1) {
        guard = find_guard_body(&z->node);
        if (guard) {
            guard_idx = find_body_index(wg, guard);
            if (guard_idx >= 0) {
                zrt->insns[ninsns].op = GEMCA_RT_GUARD_BODY;
                zrt->insns[ninsns].operand = guard_idx;
                ninsns++;
            }
        }
    }

    /* Compile the CSG tree post-order into the remaining slots. */
    compile_node(&z->node, wg, diag, zrt->insns, &ninsns);

    zrt->ninsns = ninsns;

    /*
     * Validate that the generated RPN sequence will never exceed
     * OSH_GEMCA_RT_MAX_STACK at runtime.  Simulate the stack depth by
     * replaying the instruction sequence — GUARD_BODY does not touch the
     * value stack; PUSH adds one slot; binary operators consume two and
     * produce one (net -1).  This check is O(ninsns) and runs only once at
     * setup time.
     */
    {
        int sim_sp = 0;
        int max_sp = 0;
        for (i = 0; i < ninsns; ++i) {
            switch (zrt->insns[i].op) {
            case GEMCA_RT_PUSH_BODY:
            case GEMCA_RT_PUSH_VOXEL_BODY:
                sim_sp++;
                if (sim_sp > max_sp) {
                    max_sp = sim_sp;
                }
                break;
            case GEMCA_RT_UNION:
            case GEMCA_RT_INTERSECT:
            case GEMCA_RT_DIFF:
                sim_sp--;
                break;
            default:
                break;
            }
        }
        if (max_sp > OSH_GEMCA_RT_MAX_STACK) {
            OSH_DIAG_ERRORF(diag,
                            "compile_zone: zone '%s' requires RPN stack depth %d > max %d; increase "
                            "OSH_GEMCA_RT_MAX_STACK",
                            z->name,
                            max_sp,
                            OSH_GEMCA_RT_MAX_STACK);
            free(zrt->insns);
            zrt->insns = NULL;
            return OSH_EINVAL;
        }
    }
    zrt->material_idx = z->material_idx;

    /* Detect voxel bodies for future Jacobs dispatch. */
    zrt->voxel_body_idx = -1;
    for (i = 0; i < zrt->ninsns; i++) {
        if (zrt->insns[i].op == GEMCA_RT_PUSH_VOXEL_BODY) {
            zrt->voxel_body_idx = zrt->insns[i].operand;
            break; /* at most one voxel body per zone */
        }
    }

    return OSH_OK;
}

/**
 * @brief Recursively compile a cgnode subtree into RPN instructions (post-order).
 *
 * @details
 * Leaf nodes emit a single PUSH_BODY or PUSH_VOXEL_BODY instruction.  Interior
 * nodes recurse into left and right children first, then emit the operator.
 * This post-order visit produces a valid RPN sequence: left operand, right
 * operand, operator.
 *
 * @param[in]     node    Current cgnode to compile.
 * @param[in]     wg      Cold workspace for body-pointer-to-index resolution.
 * @param[in,out] insns   Destination instruction array (must have sufficient capacity).
 * @param[in,out] ninsns  Running count of instructions written so far.
 */
static void compile_node(struct cgnode const *node,
                         struct osh_gemca_prepared const *wg,
                         struct osh_diag_sink const *diag,
                         struct gemca_rt_insn *insns,
                         int *ninsns) {
    int body_idx;

    if (node->type == _OSH_GEMCA_CGNODE_BODY) {
        if (!node->body) {
            /* Parser can produce a NULL body pointer when a body name fails to
             * resolve.  The error has already been logged; emit a sentinel
             * operand (-1) so the RPN stack stays balanced.  Evaluators treat
             * body_idx < 0 as "always outside", making the zone unreachable. */
            OSH_DIAG_ERRORF(diag, "compile_node: leaf node has NULL body pointer — zone will be unreachable");
            insns[*ninsns].op = GEMCA_RT_PUSH_BODY;
            insns[*ninsns].operand = -1;
            (*ninsns)++;
            return;
        }
        body_idx = find_body_index(wg, node->body);
        insns[*ninsns].op = (node->body->type == OSH_GEMCA_BODY_VOX) ? GEMCA_RT_PUSH_VOXEL_BODY : GEMCA_RT_PUSH_BODY;
        insns[*ninsns].operand = body_idx;
        (*ninsns)++;
        return;
    }

    /* Interior node: left subtree, right subtree, then operator. */
    compile_node(node->left, wg, diag, insns, ninsns);
    compile_node(node->right, wg, diag, insns, ninsns);

    switch (node->op) {
    case '|':
        insns[*ninsns].op = GEMCA_RT_UNION;
        break;
    case '+':
        insns[*ninsns].op = GEMCA_RT_INTERSECT;
        break;
    case '-':
        insns[*ninsns].op = GEMCA_RT_DIFF;
        break;
    default:
        OSH_DIAG_ERRORF(diag, "compile_node(): unknown CSG operator '%c'", node->op);
        insns[*ninsns].op = GEMCA_RT_UNION; /* safe fallback */
        break;
    }
    insns[*ninsns].operand = -1;
    (*ninsns)++;
}

/**
 * @brief Count leaf (body) nodes in a cgnode subtree.
 *
 * @param[in] node  Root of the subtree.
 *
 * @returns Number of leaf nodes.
 */
static int count_leaves(struct cgnode const *node) {
    if (node->type == _OSH_GEMCA_CGNODE_BODY) {
        return 1;
    }
    return count_leaves(node->left) + count_leaves(node->right);
}

/**
 * @brief Find the leftmost guard-eligible body in a cgnode tree.
 *
 * @details
 * Follows the leftmost branch from the root, descending only through
 * intersection (+) and difference (-) operators.  The first body reached on
 * this path is a valid guard: any point outside it cannot satisfy the CSG
 * expression because at least one intersection constraint fails.
 *
 * Returns NULL if the root operator is a union (|), since no single body can
 * guard a union expression.
 *
 * @param[in] node  Root cgnode of the zone's CSG tree.
 *
 * @returns Pointer to the guard body, or NULL if no single-body guard applies.
 */
static struct body const *find_guard_body(struct cgnode const *node) {
    if (node->type == _OSH_GEMCA_CGNODE_BODY) {
        return node->body;
    }
    if (node->op == '+' || node->op == '-') {
        return find_guard_body(node->left);
    }
    return NULL; /* union root: no simple guard */
}

/**
 * @brief Find the dense index of a body pointer in the cold workspace.
 *
 * @details
 * Linear search through wg->bodies[].  This is a setup-time function called
 * only during compilation; O(n) cost is acceptable.
 *
 * @param[in] wg  Cold workspace.
 * @param[in] b   Body pointer to locate.
 *
 * @returns 0-based index, or -1 if not found.
 */
static int find_body_index(struct osh_gemca_prepared const *wg, struct body const *b) {
    size_t i;
    for (i = 0; i < wg->nbodies; i++) {
        if (wg->bodies[i] == b) {
            return (int) i;
        }
    }
    return -1;
}

/* ---- RPN evaluators ------------------------------------------------------ */

/**
 * @brief Evaluate zone membership for a ray using the RPN instruction array.
 *
 * @details
 * Processes insns[] left-to-right with a small integer stack.  GUARD_BODY
 * instructions short-circuit the entire evaluation on failure (return 0
 * immediately, stack unmodified).  PUSH instructions call in_body_rt() and
 * push the result.  Binary operators pop two entries and push one combined
 * result.
 *
 * @param[in] rt  Runtime.
 * @param[in] z   Compiled zone.
 * @param[in] r   Ray to test.
 *
 * @returns 1 if the ray is inside the zone, 0 otherwise.
 */
static inline int
eval_membership(struct osh_gemca_runtime const *rt, struct gemca_rt_zone const *z, struct ray const *r) {
    int stack[OSH_GEMCA_RT_MAX_STACK];
    int sp;
    int i;
    struct gemca_rt_insn const *insn;

    /*
     * Fast path: single-instruction zones (one PUSH_BODY after the guard
     * optimisation removed guards for leaf zones) need no stack machinery.
     * This is the common case for simple body-per-zone geometry.
     */
    if (z->ninsns == 1) {
        return in_body_rt(rt, z->insns[0].operand, r);
    }

    sp = 0;

    for (i = 0; i < z->ninsns; i++) {
        insn = &z->insns[i];

        switch (insn->op) {

        case GEMCA_RT_GUARD_BODY:
            /* Fast-reject: if outside guard body, skip the entire zone. */
            if (!in_body_rt(rt, insn->operand, r)) {
                return 0;
            }
            break;

        case GEMCA_RT_PUSH_BODY:
        case GEMCA_RT_PUSH_VOXEL_BODY:
            if (sp >= OSH_GEMCA_RT_MAX_STACK) {
                return 0;
            }
            stack[sp++] = in_body_rt(rt, insn->operand, r);
            break;

        case GEMCA_RT_UNION:
            if (sp < 2) {
                return 0;
            }
            stack[sp - 2] = stack[sp - 2] || stack[sp - 1];
            sp--;
            break;

        case GEMCA_RT_INTERSECT:
            if (sp < 2) {
                return 0;
            }
            stack[sp - 2] = stack[sp - 2] && stack[sp - 1];
            sp--;
            break;

        case GEMCA_RT_DIFF:
            if (sp < 2) {
                return 0;
            }
            /* left && !right  (left = sp-2, right = sp-1 in RPN order) */
            stack[sp - 2] = stack[sp - 2] && !stack[sp - 1];
            sp--;
            break;

        default:
            return 0;
        }
    }

    return (sp > 0) ? stack[0] : 0;
}

static void eval_membership_batch_active(struct osh_gemca_runtime const *rt,
                                         struct gemca_rt_zone const *z,
                                         double const *x,
                                         double const *y,
                                         double const *zpos,
                                         double const *ux,
                                         double const *uy,
                                         double const *uz,
                                         size_t const *candidate_idx,
                                         size_t n_candidates,
                                         int *inside_out) {
    int stack[OSH_GEMCA_RT_MAX_STACK][OSH_GEMCA_RT_BODY_BATCH_CHUNK];
    int body_inside[OSH_GEMCA_RT_BODY_BATCH_CHUNK];
    size_t active_idx[OSH_GEMCA_RT_BODY_BATCH_CHUNK];
    size_t active_n;
    size_t lane;
    int sp;
    int i;
    struct gemca_rt_insn const *insn;

    if (!rt || !z || !x || !y || !zpos || !ux || !uy || !uz || !candidate_idx || !inside_out
        || n_candidates > (size_t) OSH_GEMCA_RT_BODY_BATCH_CHUNK) {
        return;
    }

    for (lane = 0; lane < n_candidates; ++lane) {
        inside_out[candidate_idx[lane]] = 0;
        active_idx[lane] = candidate_idx[lane];
    }
    active_n = n_candidates;

    if (z->ninsns == 1) {
        check_body_batch_indexed_rt(
            rt, (size_t) z->insns[0].operand, x, y, zpos, ux, uy, uz, active_idx, active_n, body_inside);
        for (lane = 0; lane < active_n; ++lane) {
            inside_out[active_idx[lane]] = body_inside[lane];
        }
        return;
    }

    sp = 0;

    for (i = 0; i < z->ninsns; ++i) {
        size_t write;

        insn = &z->insns[i];

        switch (insn->op) {
        case GEMCA_RT_GUARD_BODY:
            check_body_batch_indexed_rt(
                rt, (size_t) insn->operand, x, y, zpos, ux, uy, uz, active_idx, active_n, body_inside);
            write = 0u;
            for (lane = 0; lane < active_n; ++lane) {
                if (body_inside[lane]) {
                    active_idx[write++] = active_idx[lane];
                }
            }
            active_n = write;
            if (active_n == 0u) {
                return;
            }
            break;

        case GEMCA_RT_PUSH_BODY:
        case GEMCA_RT_PUSH_VOXEL_BODY:
            if (sp >= OSH_GEMCA_RT_MAX_STACK) {
                return;
            }
            check_body_batch_indexed_rt(
                rt, (size_t) insn->operand, x, y, zpos, ux, uy, uz, active_idx, active_n, stack[sp]);
            sp++;
            break;

        case GEMCA_RT_UNION:
            if (sp < 2) {
                return;
            }
            for (lane = 0; lane < active_n; ++lane) {
                stack[sp - 2][lane] = stack[sp - 2][lane] || stack[sp - 1][lane];
            }
            sp--;
            break;

        case GEMCA_RT_INTERSECT:
            if (sp < 2) {
                return;
            }
            for (lane = 0; lane < active_n; ++lane) {
                stack[sp - 2][lane] = stack[sp - 2][lane] && stack[sp - 1][lane];
            }
            sp--;
            break;

        case GEMCA_RT_DIFF:
            if (sp < 2) {
                return;
            }
            for (lane = 0; lane < active_n; ++lane) {
                stack[sp - 2][lane] = stack[sp - 2][lane] && !stack[sp - 1][lane];
            }
            sp--;
            break;

        default:
            return;
        }
    }

    if (sp < 1) {
        return;
    }

    for (lane = 0; lane < active_n; ++lane) {
        inside_out[active_idx[lane]] = stack[0][lane];
    }
}

/**
 * @brief Evaluate the distance to the next zone boundary at the current ray position.
 *
 * @details
 * Processes the same instruction array as eval_membership() but uses a stack
 * of (dist, is_inside) frames.  Each PUSH instruction computes both the
 * minimum positive distance to the body boundary and the inside/outside flag
 * for the current position.  Binary operators combine frames following Roth's
 * CSG ray-casting rules (Table 3, "Ray Casting for Modeling Solids", 1982):
 * all three operators use minpos(d1, d2) for the combined distance.
 *
 * GUARD_BODY instructions are no-ops here: zone membership is already confirmed
 * before get_distance is called, so every guard would trivially pass.
 *
 * @param[in]  rt         Runtime.
 * @param[in]  z          Compiled zone.
 * @param[in]  r          Ray (direction must be normalised by caller).
 * @param[out] is_inside  Set to 1 if ray is currently inside the zone, else 0.
 *
 * @returns Minimum positive distance to the nearest zone boundary.
 *
 * @todo AVX2 batch path: eval_distance has no SIMD implementation.  A future
 *       osh_gemca_runtime_get_distance_batch_avx2() would vectorise the inner
 *       body-distance loop (one body, 4 particles at a time) using _mm256_min_pd
 *       for Roth's minpos() and _mm256_blendv_pd for the CSG combine step.
 *       The step-loop in get_distance_batch must then be restructured to advance
 *       all 4 rays simultaneously.  See runtime/README.md for context.
 */
static inline double
eval_distance(struct osh_gemca_runtime const *rt, struct gemca_rt_zone const *z, struct ray const *r, int *is_inside) {
    struct dist_frame stack[OSH_GEMCA_RT_MAX_STACK];
    struct dist_frame a;
    struct dist_frame b;
    int sp;
    int i;
    struct gemca_rt_insn const *insn;

    sp = 0;

    for (i = 0; i < z->ninsns; i++) {
        insn = &z->insns[i];

        switch (insn->op) {

        case GEMCA_RT_GUARD_BODY:
            /* No-op in distance evaluation: membership already confirmed. */
            break;

        case GEMCA_RT_PUSH_BODY:
            if (sp >= OSH_GEMCA_RT_MAX_STACK) {
                *is_inside = 0;
                return 0.0;
            }
            stack[sp].is_inside = in_body_rt(rt, insn->operand, r);
            stack[sp].dist = dist_body_rt(rt, insn->operand, r);
            sp++;
            break;

        case GEMCA_RT_PUSH_VOXEL_BODY:
            if (sp >= OSH_GEMCA_RT_MAX_STACK) {
                *is_inside = 0;
                return 0.0;
            }
            stack[sp].is_inside = in_body_rt(rt, insn->operand, r);
            /* segs=NULL: CSG evaluator only needs the distance scalar.
             * Transport calls dist_voxel_body_rt() directly with its own
             * stack buffer (OSH_GEMCA_VOXEL_SEGS_STACK) to get per-voxel
             * segments for energy loss and scoring (M4/M5). */
            stack[sp].dist = dist_voxel_body_rt(rt, insn->operand, r, NULL, 0, NULL, NULL);
            sp++;
            break;

        case GEMCA_RT_UNION:
            if (sp < 2) {
                *is_inside = 0;
                return 0.0;
            }
            b = stack[--sp];
            a = stack[sp - 1];
            stack[sp - 1].is_inside = a.is_inside || b.is_inside;
            stack[sp - 1].dist = _minpos(a.dist, b.dist);
            break;

        case GEMCA_RT_INTERSECT:
            if (sp < 2) {
                *is_inside = 0;
                return 0.0;
            }
            b = stack[--sp];
            a = stack[sp - 1];
            stack[sp - 1].is_inside = a.is_inside && b.is_inside;
            stack[sp - 1].dist = _minpos(a.dist, b.dist);
            break;

        case GEMCA_RT_DIFF:
            if (sp < 2) {
                *is_inside = 0;
                return 0.0;
            }
            b = stack[--sp];
            a = stack[sp - 1];
            stack[sp - 1].is_inside = a.is_inside && !b.is_inside;
            stack[sp - 1].dist = _minpos(a.dist, b.dist);
            break;

        default:
            *is_inside = 0;
            return OSH_GEMCA_INFINITY;
        }
    }

    if (sp < 1) {
        *is_inside = 0;
        return OSH_GEMCA_INFINITY;
    }

    *is_inside = stack[0].is_inside;
    return stack[0].dist;
}

/* ---- Body evaluators ----------------------------------------------------- */

/**
 * @brief Check if a ray is inside a flat runtime body.
 *
 * @details
 * Transforms the ray to the body's local coordinate system, then checks each
 * surface in the body's flat surface slice.  Returns 0 as soon as any surface
 * reports the ray is outside (short-circuit AND).
 *
 * @param[in] rt        Runtime.
 * @param[in] body_idx  Index into rt->bodies[].
 * @param[in] r         Ray in OSH_COORD_UNIVERSE.
 *
 * @returns 1 if inside all surfaces, 0 otherwise.
 */
static inline int in_body_rt(struct osh_gemca_runtime const *rt, int body_idx, struct ray const *r) {
    struct gemca_rt_body const *b;
    struct ray tr;
    int i;

    if (body_idx < 0 || (size_t) body_idx >= rt->nbodies) {
        return 0;
    }

    b = &rt->bodies[body_idx];

    if (transform_to_local_rt(b, r, &tr) != OSH_OK) {
        return 0;
    }

    for (i = 0; i < b->nsurfs; i++) {
        if (!_check_surface_rt(&rt->surfaces[b->surf_begin + (size_t) i], &tr)) {
            return 0;
        }
    }
    return 1;
}

/**
 * @brief Compute the minimum positive distance from a ray to any surface of a body.
 *
 * @details
 * Transforms the ray to body-local coordinates, then evaluates each surface in
 * the flat slice.  Returns the smallest positive distance found, or
 * OSH_GEMCA_INFINITY if no surface is hit in the forward direction.
 *
 * @param[in] rt        Runtime.
 * @param[in] body_idx  Index into rt->bodies[].
 * @param[in] r         Ray in OSH_COORD_UNIVERSE.
 *
 * @returns Minimum positive distance to a body surface, or OSH_GEMCA_INFINITY.
 */
static inline double dist_body_rt(struct osh_gemca_runtime const *rt, int body_idx, struct ray const *r) {
    struct gemca_rt_body const *b;
    struct ray tr;
    double d;
    double min_d;
    int i;

    b = &rt->bodies[body_idx];

    if (transform_to_local_rt(b, r, &tr) != OSH_OK) {
        return OSH_GEMCA_INFINITY;
    }

    min_d = OSH_GEMCA_INFINITY;
    for (i = 0; i < b->nsurfs; i++) {
        d = _dist_surface_rt(&rt->surfaces[b->surf_begin + (size_t) i], &tr);
        if (d > 0.0 && d < min_d) {
            min_d = d;
        }
    }
    return min_d;
}

/* ---- Ray transform ------------------------------------------------------- */

/**
 * @brief Transform a ray from OSH_COORD_UNIVERSE to a body's local coordinates.
 *
 * @details
 * Mirrors _transform_to_local() in osh_gemca2_calc_zone.c and
 * osh_gemca2_dist.c but operates on a flat @ref gemca_rt_body rather than the
 * cold @ref body struct.  Three coordinate systems are handled:
 *
 *   OSH_COORD_UNIVERSE — identity; ray is copied unchanged.
 *   OSH_COORD_BCALIGN  — pure translation via the last column of t[].
 *   OSH_COORD_BZALIGN  — full affine transform via osh_ray_transform().
 *
 * @param[in]  b   Flat body (provides t[] and coord).
 * @param[in]  r   Input ray in OSH_COORD_UNIVERSE.
 * @param[out] tr  Transformed ray in body-local coordinates.
 *
 * @returns OSH_OK on success, OSH_ENOTSUP for an unknown coordinate system.
 */
static inline enum osh_status
transform_to_local_rt(struct gemca_rt_body const *b, struct ray const *r, struct ray *tr) {
    int i;
    int j;

    for (i = 0; i < 3; i++) {
        tr->p[i] = r->p[i];
        tr->cp[i] = r->cp[i];
    }
    tr->system = (int) b->coord;

    switch (b->coord) {
    case OSH_COORD_UNIVERSE:
        break;

    case OSH_COORD_BCALIGN:
        for (i = 0; i < 3; i++) {
            j = i * 4;
            tr->p[i] = r->p[i] + b->t[j + 3];
            tr->cp[i] = r->cp[i];
        }
        break;

    case OSH_COORD_BZALIGN:
        osh_ray_transform(r, tr, b->t);
        break;

    default:
        return OSH_ENOTSUP;
    }

    return OSH_OK;
}

/**
 * @brief Transform a batch of rays from OSH_COORD_UNIVERSE to one body's local coordinates.
 *
 * @details
 * Batched SoA counterpart of transform_to_local_rt().  The formulas are
 * identical to the scalar path; only the storage layout differs.
 */
static inline enum osh_status transform_to_local_batch_rt(struct gemca_rt_body const *b,
                                                          double const *x,
                                                          double const *y,
                                                          double const *z,
                                                          double const *ux,
                                                          double const *uy,
                                                          double const *uz,
                                                          size_t n,
                                                          double *tx,
                                                          double *ty,
                                                          double *tz,
                                                          double *tux,
                                                          double *tuy,
                                                          double *tuz) {
    size_t i;

    switch (b->coord) {
    case OSH_COORD_UNIVERSE:
        for (i = 0; i < n; ++i) {
            tx[i] = x[i];
            ty[i] = y[i];
            tz[i] = z[i];
            tux[i] = ux[i];
            tuy[i] = uy[i];
            tuz[i] = uz[i];
        }
        break;

    case OSH_COORD_BCALIGN:
        for (i = 0; i < n; ++i) {
            tx[i] = x[i] + b->t[3];
            ty[i] = y[i] + b->t[7];
            tz[i] = z[i] + b->t[11];
            tux[i] = ux[i];
            tuy[i] = uy[i];
            tuz[i] = uz[i];
        }
        break;

    case OSH_COORD_BZALIGN:
        for (i = 0; i < n; ++i) {
            tx[i] = x[i] * b->t[0] + y[i] * b->t[1] + z[i] * b->t[2] - b->t[3];
            ty[i] = x[i] * b->t[4] + y[i] * b->t[5] + z[i] * b->t[6] - b->t[7];
            tz[i] = x[i] * b->t[8] + y[i] * b->t[9] + z[i] * b->t[10] - b->t[11];
            tux[i] = ux[i] * b->t[0] + uy[i] * b->t[1] + uz[i] * b->t[2];
            tuy[i] = ux[i] * b->t[4] + uy[i] * b->t[5] + uz[i] * b->t[6];
            tuz[i] = ux[i] * b->t[8] + uy[i] * b->t[9] + uz[i] * b->t[10];
        }
        break;

    default:
        return OSH_ENOTSUP;
    }

    return OSH_OK;
}

/* ---- Surface evaluators -------------------------------------------------- */

/**
 * @brief Check whether one ray is on the inside of a flat surface.
 *
 * @details
 * Scalar primitive predicate shared by both the legacy scalar ray wrapper and
 * the new SoA batch surface helper.  Inputs are already in body-local
 * coordinates; no transform is applied here.
 */
static inline int _check_surface_components_rt(
    struct gemca_rt_surface const *sf, double px, double py, double pz, double ux, double uy, double uz) {
    double d;
    double dot;
    int i;
    double p[3];
    double u[3];

    p[0] = px;
    p[1] = py;
    p[2] = pz;
    u[0] = ux;
    u[1] = uy;
    u[2] = uz;

    switch (sf->type) {

    case OSH_GEMCA_SURF_SPHERE:
        d = (px * px) + (py * py) + (pz * pz) - sf->p[0];
        if (d > OSH_GEMCA_SMALL) {
            return 0;
        } else if (d < -OSH_GEMCA_SMALL) {
            return 1;
        }
        return ((px * ux) + (py * uy) + (pz * uz) < 0.0) ? 1 : 0;

    case OSH_GEMCA_SURF_ELLIPSOID:
        d = 0.0;
        for (i = 0; i < 3; i++) {
            d += (p[i] * p[i]) / sf->p[i];
        }
        if (d > 1.0 + OSH_GEMCA_SMALL) {
            return 0;
        } else if (d < 1.0 - OSH_GEMCA_SMALL) {
            return 1;
        }
        dot = 0.0;
        for (i = 0; i < 3; i++) {
            dot += (p[i] / sf->p[i]) * u[i];
        }
        return (dot < 0.0) ? 1 : 0;

    case OSH_GEMCA_SURF_CYLZ:
        d = (px * px) + (py * py) - sf->p[0];
        if (d > OSH_GEMCA_SMALL) {
            return 0;
        } else if (d < -OSH_GEMCA_SMALL) {
            return 1;
        }
        dot = (px * ux) + (py * uy);
        return (dot < 0.0) ? 1 : 0;

    case OSH_GEMCA_SURF_ELLZ:
        d = (px * px / sf->p[0]) + (py * py / sf->p[1]) - 1.0;
        if (d > OSH_GEMCA_SMALL) {
            return 0;
        } else if (d < -OSH_GEMCA_SMALL) {
            return 1;
        }
        dot = (px / sf->p[0]) * ux + (py / sf->p[1]) * uy;
        return (dot < 0.0) ? 1 : 0;

    case OSH_GEMCA_SURF_CONE:
        d = (px * px) + (py * py) - sf->p[1] * (pz * pz);
        if (d > OSH_GEMCA_SMALL) {
            return 0;
        } else if (d < -OSH_GEMCA_SMALL) {
            return 1;
        }
        dot = (px * ux) + (py * uy) - (sf->p[1] * pz * uz);
        return (dot < 0.0) ? 1 : 0;

    case OSH_GEMCA_SURF_PLANEX:
        d = sf->p[0] * px + sf->p[1];
        if (d > OSH_GEMCA_SMALL) {
            return 0;
        } else if (d < -OSH_GEMCA_SMALL) {
            return 1;
        }
        return (sf->p[0] * ux > 0.0) ? 0 : 1;

    case OSH_GEMCA_SURF_PLANEY:
        d = sf->p[0] * py + sf->p[1];
        if (d > OSH_GEMCA_SMALL) {
            return 0;
        } else if (d < -OSH_GEMCA_SMALL) {
            return 1;
        }
        return (sf->p[0] * uy > 0.0) ? 0 : 1;

    case OSH_GEMCA_SURF_PLANEZ:
        d = sf->p[0] * pz + sf->p[1];
        if (d > OSH_GEMCA_SMALL) {
            return 0;
        } else if (d < -OSH_GEMCA_SMALL) {
            return 1;
        }
        return (sf->p[0] * uz > 0.0) ? 0 : 1;

    case OSH_GEMCA_SURF_PLANE:
        d = (sf->p[0] * px) + (sf->p[1] * py) + (sf->p[2] * pz) + sf->p[3];
        if (d > OSH_GEMCA_SMALL) {
            return 0;
        } else if (d < -OSH_GEMCA_SMALL) {
            return 1;
        }
        return ((sf->p[0] * ux) + (sf->p[1] * uy) + (sf->p[2] * uz) > 0.0) ? 0 : 1;

    default:
        return 1;
    }
}

/**
 * @brief Check whether a ray is on the inside of a flat surface.
 *
 * @details
 * Dispatches to the per-type inside check.  All checks follow the same sign
 * convention as osh_gemca2_check_surface_side() in osh_gemca2_calc_surface.c:
 * a negative value of the surface implicit function means "inside".  On-
 * surface cases are disambiguated by the ray direction.
 *
 * @param[in] sf  Flat surface (type and p[] already populated at setup).
 * @param[in] r   Ray in body-local coordinates.
 *
 * @returns 1 if inside (or on boundary travelling in), 0 if outside.
 */
static inline int _check_surface_rt(struct gemca_rt_surface const *sf, struct ray const *r) {
    return _check_surface_components_rt(sf, r->p[0], r->p[1], r->p[2], r->cp[0], r->cp[1], r->cp[2]);
}

/**
 * @brief Compute the signed distance from a ray to a flat surface.
 *
 * @details
 * Dispatches to per-type distance formulas.  Returns the smallest positive
 * distance (forward intersection), 0.0 for a tangential or on-surface ray,
 * or OSH_GEMCA_INFINITY if no forward intersection exists.
 *
 * @param[in] sf  Flat surface.
 * @param[in] r   Ray in body-local coordinates (direction must be normalised).
 *
 * @returns Distance to surface along the ray direction.
 */
static inline double _dist_surface_rt(struct gemca_rt_surface const *sf, struct ray const *r) {
    switch (sf->type) {
    case OSH_GEMCA_SURF_SPHERE:
        return _dist_sphere_rt(sf->p[0], r);
    case OSH_GEMCA_SURF_ELLIPSOID:
        return _dist_ellipsoid_rt(sf->p[0], sf->p[1], sf->p[2], r);
    case OSH_GEMCA_SURF_CYLZ:
        return _dist_cyl_rt(sf->p[0], r);
    case OSH_GEMCA_SURF_ELLZ:
        return _dist_elipcyl_rt(sf->p[0], sf->p[1], r);
    case OSH_GEMCA_SURF_CONE:
        return _dist_cone_rt(sf->p[0], sf->p[1], r);
    case OSH_GEMCA_SURF_PLANEX:
        return _dist_plane_xyz_rt(0, sf, r);
    case OSH_GEMCA_SURF_PLANEY:
        return _dist_plane_xyz_rt(1, sf, r);
    case OSH_GEMCA_SURF_PLANEZ:
        return _dist_plane_xyz_rt(2, sf, r);
    case OSH_GEMCA_SURF_PLANE:
        return _dist_plane_rt(sf, r);
    default:
        return OSH_GEMCA_INFINITY;
    }
}

/* ---- Distance math (pure functions, no state) ---------------------------- */

static inline double _dist_sphere_rt(double r2, struct ray const *r) {
    double b;
    double c;
    b = 2.0 * osh_vect_dot(r->cp, r->p);
    c = osh_vect_len2(r->p) - r2;
    return _quad_solver(1.0, b, c);
}

static inline double _dist_cyl_rt(double r2, struct ray const *r) {
    double a;
    double b;
    double c;
    a = r->cp[0] * r->cp[0] + r->cp[1] * r->cp[1];
    b = 2.0 * (r->cp[0] * r->p[0] + r->cp[1] * r->p[1]);
    c = r->p[0] * r->p[0] + r->p[1] * r->p[1] - r2;
    return _quad_solver(a, b, c);
}

static inline double _dist_elipcyl_rt(double ra2, double rb2, struct ray const *r) {
    double a;
    double b;
    double c;
    a = (r->cp[0] * r->cp[0]) / ra2 + (r->cp[1] * r->cp[1]) / rb2;
    b = (r->cp[0] * r->p[0]) / ra2 + (r->cp[1] * r->p[1]) / rb2;
    c = (r->p[0] * r->p[0]) / ra2 + (r->p[1] * r->p[1]) / rb2 - 1.0;
    return _quad_solver(a, b, c);
}

static inline double _dist_cone_rt(double ra2, double rb2, struct ray const *r) {
    double a;
    double b;
    double c;
    double t;
    t = (r->p[2] - ra2) / rb2;
    a = (r->cp[0] * r->cp[0]) + (r->cp[1] * r->cp[1]) + (r->cp[2] * r->cp[2]) / rb2;
    b = (r->cp[0] * r->p[0]) + (r->cp[1] * r->p[1]) - t * r->cp[2];
    c = (r->p[0] * r->p[0]) + (r->p[1] * r->p[1]) - t * t * rb2;
    return _quad_solver(a, b, c);
}

static inline double _dist_ellipsoid_rt(double ra2, double rb2, double rc2, struct ray const *r) {
    double a;
    double b;
    double c;
    a = (r->cp[0] * r->cp[0]) / ra2 + (r->cp[1] * r->cp[1]) / rb2 + (r->cp[2] * r->cp[2]) / rc2;
    b = (r->cp[0] * r->p[0]) / ra2 + (r->cp[1] * r->p[1]) / rb2 + (r->cp[2] * r->p[2]) / rc2;
    c = (r->p[0] * r->p[0]) / ra2 + (r->p[1] * r->p[1]) / rb2 + (r->p[2] * r->p[2]) / rc2 - 1.0;
    return _quad_solver(a, b, c);
}

static inline double _dist_plane_xyz_rt(int axis, struct gemca_rt_surface const *sf, struct ray const *r) {
    double rp;
    double rcp;
    double denom;
    double d;
    rp = r->p[axis];
    rcp = r->cp[axis];
    denom = sf->p[0] * rcp;
    if (fabs(denom) < OSH_GEMCA_SMALL) {
        d = sf->p[0] * rp + sf->p[1];
        if (fabs(d) < OSH_GEMCA_SMALL) {
            return 0.0;
        }
        return OSH_GEMCA_INFINITY;
    }
    return -((sf->p[0] * rp) + sf->p[1]) / denom;
}

static inline double _dist_plane_rt(struct gemca_rt_surface const *sf, struct ray const *r) {
    double dot_ln;
    double dot_pn;
    double const *n;
    n = sf->p; /* normal: p[0..2] = A,B,C */
    dot_ln = osh_vect_dot(r->cp, n);
    if (fabs(dot_ln) < OSH_GEMCA_SMALL) {
        dot_pn = osh_vect_dot(r->p, n) + sf->p[3];
        if (fabs(dot_pn) < OSH_GEMCA_SMALL) {
            return 0.0;
        }
        return OSH_GEMCA_INFINITY;
    }
    dot_pn = osh_vect_dot(r->p, n) + sf->p[3];
    return -dot_pn / dot_ln;
}

/**
 * @brief Return the smallest positive root of a*x^2 + b*x + c = 0.
 *
 * @details
 * Returns OSH_GEMCA_INFINITY if there are no real roots or both roots are
 * non-positive.  The tangential case (discriminant == 0) produces a repeated
 * root at -b/(2a); if that root is positive it is returned as the hit distance
 * (grazing contact counts as a boundary crossing).
 *
 * @param[in] a,b,c  Quadratic coefficients.
 *
 * @returns Smallest positive root, or OSH_GEMCA_INFINITY.
 */
static inline double _quad_solver(double a, double b, double c) {
    double disc;
    double sq;
    double r1;
    double r2;

    if (fabs(a) < OSH_GEMCA_SMALL) {
        if (fabs(b) > OSH_GEMCA_SMALL) {
            r1 = -c / b;
            return (r1 > 0.0) ? r1 : OSH_GEMCA_INFINITY;
        }
        return OSH_GEMCA_INFINITY;
    }

    disc = b * b - 4.0 * a * c;
    if (disc < 0.0) {
        return OSH_GEMCA_INFINITY;
    }

    sq = sqrt(disc);
    r1 = (-b + sq) / (2.0 * a);
    r2 = (-b - sq) / (2.0 * a);
    return _minpos(r1, r2);
}

/**
 * @brief Return the smallest positive value of two doubles, or 0 if both are non-positive.
 *
 * @param[in] a,b  Two values to compare.
 *
 * @returns Smallest positive value, or 0.0.
 */
static inline double _minpos(double a, double b) {
    if (a > 0.0) {
        if (b > 0.0) {
            return _RT_MIN(a, b);
        }
        return a;
    }
    if (b > 0.0) {
        return b;
    }
    return 0.0;
}
