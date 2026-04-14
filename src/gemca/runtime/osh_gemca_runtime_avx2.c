/**
 * @file osh_gemca_runtime_avx2.c
 * @brief AVX2+FMA batch zone-membership implementation.
 *
 * @details
 * This translation unit is compiled with -mavx2 -mfma and must not be
 * included in a build where those flags are absent.  The public entry point
 * osh_gemca_runtime_get_zone_batch_avx2() is called from the scalar fallback
 * in osh_gemca_runtime.c after a runtime __builtin_cpu_supports("avx2") check.
 *
 * Algorithm
 * ---------
 * Zone-outer, 4-particle-inner loop:
 *   For each group of 4 particles, load six __m256d vectors (x,y,z,ux,uy,uz).
 *   Iterate zones until all 4 particles have been resolved (zone index found or
 *   exhausted).  For each zone, _eval_membership_avx2() evaluates the RPN
 *   instruction array 4-wide and returns a __m256i lane mask (all-bits-set =
 *   inside, all-zero = outside).  Newly-resolved particles are recorded via a
 *   4-bit scalar bitmask; the zone index is written with a scalar store.
 *   A scalar tail handles the final 0-3 particles that do not form a full group.
 *
 * Surface sign convention
 * -----------------------
 * _surf_sign_avx2(vd, vdot) returns all-bits-set for "inside":
 *   - vd < -SMALL                        → inside
 *   - |vd| ≤ SMALL  AND  vdot < 0       → inside (boundary: entering surface)
 *   - vd > SMALL                         → outside
 *   - |vd| ≤ SMALL  AND  vdot ≥ 0       → outside
 * This matches the scalar _check_surface_components_rt() logic exactly, with
 * a measure-zero difference at the vdot = 0 boundary (scalar says inside,
 * AVX2 says outside — irrelevant in practice).
 *
 * GPU migration note
 * ------------------
 * The structure of _eval_membership_avx2() maps almost directly to a GPU
 * kernel body: replace __m256d with scalar double (the GPU does its own SIMD
 * across threads in a warp), replace __m256i boolean lanes with a single int,
 * and launch one thread per particle.  The RPN stack fits in registers at
 * OSH_GEMCA_RT_MAX_STACK = 32 slots.  The only prerequisite is flattening
 * zones[j].insns into a device-accessible insns_flat[] array — see
 * runtime/README.md.
 */

#include <immintrin.h>
#include <stddef.h>

#include "common/osh_coord.h"
#include "common/osh_ray.h"
#include "gemca/osh_gemca2.h"
#include "gemca/osh_gemca2_defines.h"
#include "gemca/runtime/osh_gemca_runtime.h"

/* ---- Sign helper ---------------------------------------------------------- */

/**
 * @brief Convert (implicit-function value, dot-product) to a 4-lane inside mask.
 *
 * @param vd    Implicit function value for 4 particles.
 * @param vdot  Direction dot product (outward normal · direction) for 4 particles.
 * @returns     __m256i: all-bits-set lanes = inside, all-zero = outside.
 */
static inline __m256i _surf_sign_avx2(__m256d vd, __m256d vdot) {
    __m256d const vsmall = _mm256_set1_pd(OSH_GEMCA_SMALL);
    __m256d const vnsmall = _mm256_set1_pd(-OSH_GEMCA_SMALL);
    __m256d const vzero = _mm256_setzero_pd();

    /* d < -SMALL → clearly inside */
    __m256d m_inside = _mm256_cmp_pd(vd, vnsmall, _CMP_LT_OQ);

    /* |d| ≤ SMALL: use direction to break tie */
    __m256d m_le = _mm256_cmp_pd(vd, vsmall, _CMP_LE_OQ);
    __m256d m_ge = _mm256_cmp_pd(vd, vnsmall, _CMP_GE_OQ);
    __m256d m_on = _mm256_and_pd(m_le, m_ge);                  /* |d| ≤ SMALL   */
    __m256d m_dot_lt = _mm256_cmp_pd(vdot, vzero, _CMP_LT_OQ); /* dot < 0  */

    m_inside = _mm256_or_pd(m_inside, _mm256_and_pd(m_on, m_dot_lt));

    return _mm256_castpd_si256(m_inside);
}

/* ---- Surface evaluator ---------------------------------------------------- */

/**
 * @brief Evaluate zone-membership for one surface for 4 particles simultaneously.
 *
 * @param sf   Flat runtime surface (type + parameters).
 * @param tx   Transformed x-positions for 4 particles (already in body-local coords).
 * @param ty   Transformed y-positions.
 * @param tz   Transformed z-positions.
 * @param tux  Transformed x-directions.
 * @param tuy  Transformed y-directions.
 * @param tuz  Transformed z-directions.
 * @returns    __m256i lane mask: all-bits-set = inside surface, all-zero = outside.
 */
static __m256i _check_surface_avx2(
    struct gemca_rt_surface const *sf, __m256d tx, __m256d ty, __m256d tz, __m256d tux, __m256d tuy, __m256d tuz) {
    double const *p = sf->p;
    __m256d vd, vdot;

    switch (sf->type) {

    case OSH_GEMCA_SURF_SPHERE: {
        /* f = x² + y² + z² - p[0],  dot = x·ux + y·uy + z·uz */
        __m256d vp0 = _mm256_set1_pd(p[0]);
        vd = _mm256_sub_pd(_mm256_fmadd_pd(tx, tx, _mm256_fmadd_pd(ty, ty, _mm256_mul_pd(tz, tz))), vp0);
        vdot = _mm256_fmadd_pd(tx, tux, _mm256_fmadd_pd(ty, tuy, _mm256_mul_pd(tz, tuz)));
        return _surf_sign_avx2(vd, vdot);
    }

    case OSH_GEMCA_SURF_ELLIPSOID: {
        /* f = x²/p0 + y²/p1 + z²/p2 - 1,  dot = (x/p0)·ux + (y/p1)·uy + (z/p2)·uz */
        double inv0 = 1.0 / p[0];
        double inv1 = 1.0 / p[1];
        double inv2 = 1.0 / p[2];
        __m256d vi0 = _mm256_set1_pd(inv0);
        __m256d vi1 = _mm256_set1_pd(inv1);
        __m256d vi2 = _mm256_set1_pd(inv2);
        vd = _mm256_sub_pd(
            _mm256_fmadd_pd(_mm256_mul_pd(tx, tx),
                            vi0,
                            _mm256_fmadd_pd(_mm256_mul_pd(ty, ty), vi1, _mm256_mul_pd(_mm256_mul_pd(tz, tz), vi2))),
            _mm256_set1_pd(1.0));
        vdot =
            _mm256_fmadd_pd(_mm256_mul_pd(tx, vi0),
                            tux,
                            _mm256_fmadd_pd(_mm256_mul_pd(ty, vi1), tuy, _mm256_mul_pd(_mm256_mul_pd(tz, vi2), tuz)));
        return _surf_sign_avx2(vd, vdot);
    }

    case OSH_GEMCA_SURF_CYLZ: {
        /* f = x² + y² - p[0],  dot = x·ux + y·uy */
        __m256d vp0 = _mm256_set1_pd(p[0]);
        vd = _mm256_sub_pd(_mm256_fmadd_pd(tx, tx, _mm256_mul_pd(ty, ty)), vp0);
        vdot = _mm256_fmadd_pd(tx, tux, _mm256_mul_pd(ty, tuy));
        return _surf_sign_avx2(vd, vdot);
    }

    case OSH_GEMCA_SURF_ELLZ: {
        /* f = x²/p0 + y²/p1 - 1,  dot = (x/p0)·ux + (y/p1)·uy */
        double inv0 = 1.0 / p[0];
        double inv1 = 1.0 / p[1];
        __m256d vi0 = _mm256_set1_pd(inv0);
        __m256d vi1 = _mm256_set1_pd(inv1);
        vd = _mm256_sub_pd(_mm256_fmadd_pd(_mm256_mul_pd(tx, tx), vi0, _mm256_mul_pd(_mm256_mul_pd(ty, ty), vi1)),
                           _mm256_set1_pd(1.0));
        vdot = _mm256_fmadd_pd(_mm256_mul_pd(tx, vi0), tux, _mm256_mul_pd(_mm256_mul_pd(ty, vi1), tuy));
        return _surf_sign_avx2(vd, vdot);
    }

    case OSH_GEMCA_SURF_CONE: {
        /* f = x² + y² - p[1]·z²,  dot = x·ux + y·uy - p[1]·z·uz */
        __m256d vp1 = _mm256_set1_pd(p[1]);
        vd = _mm256_sub_pd(_mm256_fmadd_pd(tx, tx, _mm256_mul_pd(ty, ty)), _mm256_mul_pd(vp1, _mm256_mul_pd(tz, tz)));
        vdot =
            _mm256_sub_pd(_mm256_fmadd_pd(tx, tux, _mm256_mul_pd(ty, tuy)), _mm256_mul_pd(vp1, _mm256_mul_pd(tz, tuz)));
        return _surf_sign_avx2(vd, vdot);
    }

    case OSH_GEMCA_SURF_PLANEX: {
        /* f = p[0]·x + p[1],  dot = p[0]·ux */
        __m256d vp0 = _mm256_set1_pd(p[0]);
        __m256d vp1 = _mm256_set1_pd(p[1]);
        vd = _mm256_fmadd_pd(vp0, tx, vp1);
        vdot = _mm256_mul_pd(vp0, tux);
        return _surf_sign_avx2(vd, vdot);
    }

    case OSH_GEMCA_SURF_PLANEY: {
        /* f = p[0]·y + p[1],  dot = p[0]·uy */
        __m256d vp0 = _mm256_set1_pd(p[0]);
        __m256d vp1 = _mm256_set1_pd(p[1]);
        vd = _mm256_fmadd_pd(vp0, ty, vp1);
        vdot = _mm256_mul_pd(vp0, tuy);
        return _surf_sign_avx2(vd, vdot);
    }

    case OSH_GEMCA_SURF_PLANEZ: {
        /* f = p[0]·z + p[1],  dot = p[0]·uz */
        __m256d vp0 = _mm256_set1_pd(p[0]);
        __m256d vp1 = _mm256_set1_pd(p[1]);
        vd = _mm256_fmadd_pd(vp0, tz, vp1);
        vdot = _mm256_mul_pd(vp0, tuz);
        return _surf_sign_avx2(vd, vdot);
    }

    case OSH_GEMCA_SURF_PLANE: {
        /* f = p[0]·x + p[1]·y + p[2]·z + p[3],  dot = p[0]·ux + p[1]·uy + p[2]·uz */
        __m256d vp0 = _mm256_set1_pd(p[0]);
        __m256d vp1 = _mm256_set1_pd(p[1]);
        __m256d vp2 = _mm256_set1_pd(p[2]);
        __m256d vp3 = _mm256_set1_pd(p[3]);
        vd = _mm256_fmadd_pd(vp0, tx, _mm256_fmadd_pd(vp1, ty, _mm256_fmadd_pd(vp2, tz, vp3)));
        vdot = _mm256_fmadd_pd(vp0, tux, _mm256_fmadd_pd(vp1, tuy, _mm256_mul_pd(vp2, tuz)));
        return _surf_sign_avx2(vd, vdot);
    }

    default:
        return _mm256_setzero_si256();
    }
}

/* ---- Body evaluator ------------------------------------------------------- */

/**
 * @brief Test whether 4 particles are inside body[body_idx].
 *
 * @details
 * Transforms the 4-wide position and direction vectors into the body's local
 * coordinate system, then ANDs the inside result of every surface.  Early exit
 * if all 4 particles are already outside (testz returns non-zero).
 *
 * @param rt        Runtime.
 * @param body_idx  Index into rt->bodies[].
 * @param vx … vuz  4-wide position and direction vectors in universe coordinates.
 * @returns         __m256i lane mask: all-bits-set = inside body, all-zero = outside.
 */
static __m256i _in_body_avx2(struct gemca_runtime const *rt,
                             int body_idx,
                             __m256d vx,
                             __m256d vy,
                             __m256d vz,
                             __m256d vux,
                             __m256d vuy,
                             __m256d vuz) {
    struct gemca_rt_body const *b;
    __m256d tx, ty, tz, tux, tuy, tuz;
    __m256i inside;
    int is;

    if (body_idx < 0 || (size_t) body_idx >= rt->nbodies) {
        return _mm256_setzero_si256(); /* sentinel / unresolved body → always outside */
    }

    b = &rt->bodies[body_idx];

    switch (b->coord) {
    case OSH_COORD_UNIVERSE:
        tx = vx;
        ty = vy;
        tz = vz;
        tux = vux;
        tuy = vuy;
        tuz = vuz;
        break;

    case OSH_COORD_BCALIGN: {
        /* Translation only: p_local = p_universe + t[3,7,11] */
        __m256d dt3 = _mm256_set1_pd(b->t[3]);
        __m256d dt7 = _mm256_set1_pd(b->t[7]);
        __m256d dt11 = _mm256_set1_pd(b->t[11]);
        tx = _mm256_add_pd(vx, dt3);
        ty = _mm256_add_pd(vy, dt7);
        tz = _mm256_add_pd(vz, dt11);
        tux = vux;
        tuy = vuy;
        tuz = vuz;
        break;
    }

    case OSH_COORD_BZALIGN: {
        /*
         * Full 3×4 affine (row-major 4×4):
         *   p_local = R · p_universe − t_col
         *   tx = x·t[0] + y·t[1] + z·t[2] − t[3]
         *   ty = x·t[4] + y·t[5] + z·t[6] − t[7]
         *   tz = x·t[8] + y·t[9] + z·t[10] − t[11]
         *
         * Direction vector (no translation):
         *   tux = ux·t[0] + uy·t[1] + uz·t[2]
         *   tuy = ux·t[4] + uy·t[5] + uz·t[6]
         *   tuz = ux·t[8] + uy·t[9] + uz·t[10]
         */
        __m256d b0 = _mm256_set1_pd(b->t[0]);
        __m256d b1 = _mm256_set1_pd(b->t[1]);
        __m256d b2 = _mm256_set1_pd(b->t[2]);
        __m256d b3 = _mm256_set1_pd(b->t[3]);
        __m256d b4 = _mm256_set1_pd(b->t[4]);
        __m256d b5 = _mm256_set1_pd(b->t[5]);
        __m256d b6 = _mm256_set1_pd(b->t[6]);
        __m256d b7 = _mm256_set1_pd(b->t[7]);
        __m256d b8 = _mm256_set1_pd(b->t[8]);
        __m256d b9 = _mm256_set1_pd(b->t[9]);
        __m256d b10 = _mm256_set1_pd(b->t[10]);
        __m256d b11 = _mm256_set1_pd(b->t[11]);

        tx = _mm256_fmadd_pd(vx, b0, _mm256_fmadd_pd(vy, b1, _mm256_fmsub_pd(vz, b2, b3)));
        ty = _mm256_fmadd_pd(vx, b4, _mm256_fmadd_pd(vy, b5, _mm256_fmsub_pd(vz, b6, b7)));
        tz = _mm256_fmadd_pd(vx, b8, _mm256_fmadd_pd(vy, b9, _mm256_fmsub_pd(vz, b10, b11)));
        tux = _mm256_fmadd_pd(vux, b0, _mm256_fmadd_pd(vuy, b1, _mm256_mul_pd(vuz, b2)));
        tuy = _mm256_fmadd_pd(vux, b4, _mm256_fmadd_pd(vuy, b5, _mm256_mul_pd(vuz, b6)));
        tuz = _mm256_fmadd_pd(vux, b8, _mm256_fmadd_pd(vuy, b9, _mm256_mul_pd(vuz, b10)));
        break;
    }

    default:
        return _mm256_setzero_si256();
    }

    /* Start with all-inside, then AND each surface result */
    inside = _mm256_set1_epi64x(-1LL);

    for (is = 0; is < b->nsurfs; ++is) {
        inside = _mm256_and_si256(
            inside, _check_surface_avx2(&rt->surfaces[b->surf_begin + (size_t) is], tx, ty, tz, tux, tuy, tuz));
        /* Early exit: all 4 particles are outside this body */
        if (_mm256_testz_si256(inside, inside)) {
            return _mm256_setzero_si256();
        }
    }

    return inside;
}

/* ---- RPN zone evaluator --------------------------------------------------- */

/**
 * @brief Evaluate a zone's RPN membership expression for 4 particles simultaneously.
 *
 * @details
 * Mirrors eval_membership() in osh_gemca_runtime.c but operates on 4 particles
 * at once using a __m256i stack (each 64-bit lane: all-bits-set = true/inside,
 * all-zero = false/outside).
 *
 * GUARD_BODY: if all 4 particles fail the guard, return all-zero immediately.
 * Individual per-lane guard failures are handled implicitly — the guard body also
 * appears as a PUSH_BODY in the RPN expression and will produce a 0 for outside
 * particles, making the final zone result 0.
 *
 * @param rt       Runtime.
 * @param z        Compiled zone.
 * @param vx … vuz 4-wide position + direction in universe coordinates.
 * @returns        __m256i: all-bits-set lanes = inside zone, all-zero = outside.
 */
static __m256i _eval_membership_avx2(struct gemca_runtime const *rt,
                                     struct gemca_rt_zone const *z,
                                     __m256d vx,
                                     __m256d vy,
                                     __m256d vz,
                                     __m256d vux,
                                     __m256d vuy,
                                     __m256d vuz) {
    __m256i vstack[OSH_GEMCA_RT_MAX_STACK];
    int sp = 0;
    int i;
    struct gemca_rt_insn const *insn;

    /* Fast path for single-body zones (no guard, no operators) */
    if (z->ninsns == 1) {
        return _in_body_avx2(rt, z->insns[0].operand, vx, vy, vz, vux, vuy, vuz);
    }

    for (i = 0; i < z->ninsns; ++i) {
        insn = &z->insns[i];

        switch (insn->op) {

        case GEMCA_RT_GUARD_BODY: {
            __m256i g = _in_body_avx2(rt, insn->operand, vx, vy, vz, vux, vuy, vuz);
            /* All 4 particles outside the guard body → reject whole group */
            if (_mm256_testz_si256(g, g)) {
                return _mm256_setzero_si256();
            }
            break;
        }

        case GEMCA_RT_PUSH_BODY:
        case GEMCA_RT_PUSH_VOXEL_BODY:
            if (sp >= OSH_GEMCA_RT_MAX_STACK) {
                return _mm256_setzero_si256(); /* stack overflow → outside */
            }
            vstack[sp++] = _in_body_avx2(rt, insn->operand, vx, vy, vz, vux, vuy, vuz);
            break;

        case GEMCA_RT_UNION:
            if (sp < 2) {
                return _mm256_setzero_si256();
            }
            vstack[sp - 2] = _mm256_or_si256(vstack[sp - 2], vstack[sp - 1]);
            sp--;
            break;

        case GEMCA_RT_INTERSECT:
            if (sp < 2) {
                return _mm256_setzero_si256();
            }
            vstack[sp - 2] = _mm256_and_si256(vstack[sp - 2], vstack[sp - 1]);
            sp--;
            break;

        case GEMCA_RT_DIFF:
            if (sp < 2) {
                return _mm256_setzero_si256();
            }
            /* left AND NOT right */
            vstack[sp - 2] = _mm256_andnot_si256(vstack[sp - 1], vstack[sp - 2]);
            sp--;
            break;

        default:
            return _mm256_setzero_si256();
        }
    }

    return (sp > 0) ? vstack[0] : _mm256_setzero_si256();
}

/* ---- Public entry point --------------------------------------------------- */

/**
 * @brief AVX2+FMA batch zone lookup.
 *
 * @details
 * Processes @p n particles in groups of 4 using 256-bit wide SIMD.  Each
 * group iterates over zones until all 4 particles are resolved or all zones
 * are exhausted.  The remaining 0-3 particles are handled by a scalar loop
 * using osh_gemca_runtime_get_zone().
 *
 * This function is only compiled when the build system detects -mavx2 -mfma
 * support and defines OSH_GEMCA_RUNTIME_HAVE_AVX2=1.  It is called via a
 * runtime dispatch in osh_gemca_runtime_get_zone_batch().
 */
void osh_gemca_runtime_get_zone_batch_avx2(struct gemca_runtime const *rt,
                                           double const *x,
                                           double const *y,
                                           double const *z,
                                           double const *ux,
                                           double const *uy,
                                           double const *uz,
                                           size_t n,
                                           size_t *zone_out) {
    size_t base;
    size_t j;

    /* --- 4-wide AVX2 loop --- */
    for (base = 0; base + 4 <= n; base += 4) {
        __m256d vx = _mm256_loadu_pd(x + base);
        __m256d vy = _mm256_loadu_pd(y + base);
        __m256d vz = _mm256_loadu_pd(z + base);
        __m256d vux = _mm256_loadu_pd(ux + base);
        __m256d vuy = _mm256_loadu_pd(uy + base);
        __m256d vuz = _mm256_loadu_pd(uz + base);

        zone_out[base + 0] = OSH_GEMCA_ZONE_INDEX_INVALID;
        zone_out[base + 1] = OSH_GEMCA_ZONE_INDEX_INVALID;
        zone_out[base + 2] = OSH_GEMCA_ZONE_INDEX_INVALID;
        zone_out[base + 3] = OSH_GEMCA_ZONE_INDEX_INVALID;

        /*
         * 4-bit scalar bitmask: bit k = 1 means particle (base+k) has been
         * assigned a zone and should no longer update zone_out.
         */
        int resolved = 0;

        for (j = 0; j < rt->nzones; ++j) {
            __m256i inside;
            int m;
            int newly;

            if (resolved == 0xF) {
                break; /* All 4 particles resolved */
            }

            inside = _eval_membership_avx2(rt, &rt->zones[j], vx, vy, vz, vux, vuy, vuz);

            /*
             * _mm256_movemask_pd reads the sign bit of each 64-bit double lane.
             * Because inside masks are all-bits-set (sign bit = 1) or all-zero
             * (sign bit = 0), movemask gives a 4-bit integer with one bit per lane.
             */
            m = _mm256_movemask_pd(_mm256_castsi256_pd(inside));

            /* Newly resolved: inside this zone AND not yet assigned */
            newly = m & ~resolved;

            if (newly & 1) {
                zone_out[base + 0] = j;
            }
            if (newly & 2) {
                zone_out[base + 1] = j;
            }
            if (newly & 4) {
                zone_out[base + 2] = j;
            }
            if (newly & 8) {
                zone_out[base + 3] = j;
            }

            resolved |= newly;
        }
    }

    /* --- Scalar tail for remaining 0-3 particles --- */
    for (; base < n; ++base) {
        struct ray r;
        r.p[0] = x[base];
        r.p[1] = y[base];
        r.p[2] = z[base];
        r.cp[0] = ux[base];
        r.cp[1] = uy[base];
        r.cp[2] = uz[base];
        r.system = OSH_COORD_UNIVERSE;
        zone_out[base] = osh_gemca_runtime_get_zone(rt, &r);
    }
}
