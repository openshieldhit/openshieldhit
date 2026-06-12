#include "gemca/runtime/osh_gemca_runtime_accel.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "gemca/osh_gemca2.h"
#include "gemca/osh_gemca2_defines.h"

/* ---- Local constants ------------------------------------------------------ */

/* Target number of grid cells per zone; total is clamped below. */
#define ACCEL_CELLS_PER_ZONE 8
#define ACCEL_MIN_CELLS 64
#define ACCEL_MAX_CELLS_PER_AXIS 128
/* Hard cap on the specialized instruction pool (entries). */
#define ACCEL_MAX_POOL_INSNS (8u * 1024u * 1024u)

/* ---- AABB helpers --------------------------------------------------------- */

static void aabb_set_unbounded(struct gemca_rt_aabb *bb) {
    int i;
    for (i = 0; i < 3; i++) {
        bb->lo[i] = -HUGE_VAL;
        bb->hi[i] = HUGE_VAL;
    }
}

static int aabb_is_bounded(struct gemca_rt_aabb const *bb) {
    int i;
    for (i = 0; i < 3; i++) {
        if (!isfinite(bb->lo[i]) || !isfinite(bb->hi[i])) {
            return 0;
        }
    }
    return 1;
}

static void aabb_intersect(struct gemca_rt_aabb *a, struct gemca_rt_aabb const *b) {
    int i;
    for (i = 0; i < 3; i++) {
        if (b->lo[i] > a->lo[i]) {
            a->lo[i] = b->lo[i];
        }
        if (b->hi[i] < a->hi[i]) {
            a->hi[i] = b->hi[i];
        }
    }
}

static void aabb_union(struct gemca_rt_aabb *a, struct gemca_rt_aabb const *b) {
    int i;
    for (i = 0; i < 3; i++) {
        if (b->lo[i] < a->lo[i]) {
            a->lo[i] = b->lo[i];
        }
        if (b->hi[i] > a->hi[i]) {
            a->hi[i] = b->hi[i];
        }
    }
}

static int aabb_overlaps(struct gemca_rt_aabb const *a, struct gemca_rt_aabb const *b) {
    int i;
    for (i = 0; i < 3; i++) {
        if (a->lo[i] > b->hi[i] || a->hi[i] < b->lo[i]) {
            return 0;
        }
    }
    return 1;
}

/*
 * Conservative inflation: the runtime membership predicates accept points
 * within OSH_GEMCA_SMALL of a surface (implicit-function epsilon), and the
 * universe-frame corner transform rounds.  Pad generously; the pad is tiny
 * relative to any realistic geometry feature and culling stays effective.
 */
static void aabb_inflate(struct gemca_rt_aabb *bb) {
    int i;
    double pad;
    for (i = 0; i < 3; i++) {
        if (isfinite(bb->lo[i])) {
            pad = 1e-6 + 1e-9 * fabs(bb->lo[i]);
            bb->lo[i] -= pad;
        }
        if (isfinite(bb->hi[i])) {
            pad = 1e-6 + 1e-9 * fabs(bb->hi[i]);
            bb->hi[i] += pad;
        }
    }
}

/* Tighten one axis of bb to [lo, hi] (either may be infinite). */
static void aabb_clamp_axis(struct gemca_rt_aabb *bb, int axis, double lo, double hi) {
    if (lo > bb->lo[axis]) {
        bb->lo[axis] = lo;
    }
    if (hi < bb->hi[axis]) {
        bb->hi[axis] = hi;
    }
}

/* ---- Body AABB derivation ------------------------------------------------- */

/*
 * Accumulate the axis bounds implied by one surface (body-local frame).
 * A body is the intersection of all its surface half-spaces, so each
 * surface can only tighten the box.  Surfaces that do not constrain an
 * axis (cones, arbitrary planes) leave it untouched.
 */
static void surface_local_bounds(struct gemca_rt_surface const *sf, struct gemca_rt_aabb *bb) {
    double r;

    switch (sf->type) {
    case OSH_GEMCA_SURF_SPHERE:
        r = sqrt(sf->p[0]);
        aabb_clamp_axis(bb, 0, -r, r);
        aabb_clamp_axis(bb, 1, -r, r);
        aabb_clamp_axis(bb, 2, -r, r);
        break;

    case OSH_GEMCA_SURF_ELLIPSOID:
        aabb_clamp_axis(bb, 0, -sqrt(sf->p[0]), sqrt(sf->p[0]));
        aabb_clamp_axis(bb, 1, -sqrt(sf->p[1]), sqrt(sf->p[1]));
        aabb_clamp_axis(bb, 2, -sqrt(sf->p[2]), sqrt(sf->p[2]));
        break;

    case OSH_GEMCA_SURF_CYLZ:
        r = sqrt(sf->p[0]);
        aabb_clamp_axis(bb, 0, -r, r);
        aabb_clamp_axis(bb, 1, -r, r);
        break;

    case OSH_GEMCA_SURF_ELLZ:
        aabb_clamp_axis(bb, 0, -sqrt(sf->p[0]), sqrt(sf->p[0]));
        aabb_clamp_axis(bb, 1, -sqrt(sf->p[1]), sqrt(sf->p[1]));
        break;

    case OSH_GEMCA_SURF_PLANEX:
        /* inside: p0*x + p1 <= 0 */
        if (sf->p[0] > 0.0) {
            aabb_clamp_axis(bb, 0, -HUGE_VAL, -sf->p[1] / sf->p[0]);
        } else if (sf->p[0] < 0.0) {
            aabb_clamp_axis(bb, 0, -sf->p[1] / sf->p[0], HUGE_VAL);
        }
        break;

    case OSH_GEMCA_SURF_PLANEY:
        if (sf->p[0] > 0.0) {
            aabb_clamp_axis(bb, 1, -HUGE_VAL, -sf->p[1] / sf->p[0]);
        } else if (sf->p[0] < 0.0) {
            aabb_clamp_axis(bb, 1, -sf->p[1] / sf->p[0], HUGE_VAL);
        }
        break;

    case OSH_GEMCA_SURF_PLANEZ:
        if (sf->p[0] > 0.0) {
            aabb_clamp_axis(bb, 2, -HUGE_VAL, -sf->p[1] / sf->p[0]);
        } else if (sf->p[0] < 0.0) {
            aabb_clamp_axis(bb, 2, -sf->p[1] / sf->p[0], HUGE_VAL);
        }
        break;

    case OSH_GEMCA_SURF_CONE:
    case OSH_GEMCA_SURF_PLANE:
    default:
        /* No single-axis bound derivable in isolation. */
        break;
    }
}

/* Second pass: a cone (x^2 + y^2 <= p1*z^2) bounds x and y once z is bounded. */
static void cone_lateral_bounds(struct gemca_rt_surface const *sf, struct gemca_rt_aabb *bb) {
    double zmax;
    double r;

    if (sf->type != OSH_GEMCA_SURF_CONE) {
        return;
    }
    if (!isfinite(bb->lo[2]) || !isfinite(bb->hi[2])) {
        return;
    }
    zmax = fabs(bb->lo[2]) > fabs(bb->hi[2]) ? fabs(bb->lo[2]) : fabs(bb->hi[2]);
    r = sqrt(sf->p[1]) * zmax;
    aabb_clamp_axis(bb, 0, -r, r);
    aabb_clamp_axis(bb, 1, -r, r);
}

/* Check that the 3x3 rotation part of a BZALIGN transform is orthonormal,
 * so its transpose is an exact inverse for mapping the local box back to
 * the universe frame. */
static int rotation_is_orthonormal(double const *t) {
    double dot;
    int row_a;
    int row_b;
    int k;

    for (row_a = 0; row_a < 3; row_a++) {
        for (row_b = 0; row_b < 3; row_b++) {
            dot = 0.0;
            for (k = 0; k < 3; k++) {
                dot += t[row_a * 4 + k] * t[row_b * 4 + k];
            }
            if (fabs(dot - ((row_a == row_b) ? 1.0 : 0.0)) > 1e-9) {
                return 0;
            }
        }
    }
    return 1;
}

/*
 * Compute the universe-frame AABB of one body.  The local box is derived
 * from the body's surface intersection, then mapped through the inverse of
 * the body transform.  Any axis that stays unbounded (or a transform we
 * cannot invert safely) marks the body unbounded — it is then never culled.
 */
static void
body_universe_aabb(struct osh_gemca_runtime const *rt, size_t body_idx, struct gemca_rt_aabb *bb, uint8_t *bounded) {
    struct gemca_rt_body const *b = &rt->bodies[body_idx];
    struct gemca_rt_aabb local;
    double corner_l[3];
    double corner_u[3];
    double tv[3];
    int i;
    int cx;
    int cy;
    int cz;

    aabb_set_unbounded(&local);
    for (i = 0; i < b->nsurfs; i++) {
        surface_local_bounds(&rt->surfaces[b->surf_begin + (size_t) i], &local);
    }
    for (i = 0; i < b->nsurfs; i++) {
        cone_lateral_bounds(&rt->surfaces[b->surf_begin + (size_t) i], &local);
    }

    switch (b->coord) {
    case OSH_COORD_UNIVERSE:
        *bb = local;
        break;

    case OSH_COORD_BCALIGN:
        /* local = universe + (t3, t7, t11) => universe = local - translation */
        *bb = local;
        for (i = 0; i < 3; i++) {
            bb->lo[i] -= b->t[i * 4 + 3];
            bb->hi[i] -= b->t[i * 4 + 3];
        }
        break;

    case OSH_COORD_BZALIGN:
        /* local = R*universe - tv => universe = R^T * (local + tv) */
        if (!aabb_is_bounded(&local) || !rotation_is_orthonormal(b->t)) {
            aabb_set_unbounded(bb);
            break;
        }
        tv[0] = b->t[3];
        tv[1] = b->t[7];
        tv[2] = b->t[11];
        aabb_set_unbounded(bb);
        bb->lo[0] = HUGE_VAL;
        bb->lo[1] = HUGE_VAL;
        bb->lo[2] = HUGE_VAL;
        bb->hi[0] = -HUGE_VAL;
        bb->hi[1] = -HUGE_VAL;
        bb->hi[2] = -HUGE_VAL;
        for (cx = 0; cx < 2; cx++) {
            for (cy = 0; cy < 2; cy++) {
                for (cz = 0; cz < 2; cz++) {
                    corner_l[0] = (cx ? local.hi[0] : local.lo[0]) + tv[0];
                    corner_l[1] = (cy ? local.hi[1] : local.lo[1]) + tv[1];
                    corner_l[2] = (cz ? local.hi[2] : local.lo[2]) + tv[2];
                    for (i = 0; i < 3; i++) {
                        /* column i of R == row i of R^T */
                        corner_u[i] = corner_l[0] * b->t[0 * 4 + i] + corner_l[1] * b->t[1 * 4 + i]
                                      + corner_l[2] * b->t[2 * 4 + i];
                        if (corner_u[i] < bb->lo[i]) {
                            bb->lo[i] = corner_u[i];
                        }
                        if (corner_u[i] > bb->hi[i]) {
                            bb->hi[i] = corner_u[i];
                        }
                    }
                }
            }
        }
        break;

    default:
        aabb_set_unbounded(bb);
        break;
    }

    aabb_inflate(bb);
    *bounded = (uint8_t) aabb_is_bounded(bb);
}

/* ---- Zone AABB ------------------------------------------------------------- */

/*
 * Fold the zone's RPN program over AABBs: PUSH -> body box, UNION -> box
 * union, INTERSECT -> box intersection, DIFF -> left box (a difference is a
 * subset of its left operand).  A trailing intersection with the guard body's
 * box is valid because the guard is an intersection term of the expression.
 */
static int zone_aabb(struct osh_gemca_runtime const *rt,
                     struct osh_gemca_accel const *accel,
                     struct gemca_rt_zone const *z,
                     struct gemca_rt_aabb *out) {
    struct gemca_rt_aabb stack[OSH_GEMCA_RT_MAX_STACK];
    int sp = 0;
    int i;
    int guard_idx = -1;
    struct gemca_rt_insn const *insn;

    for (i = 0; i < z->ninsns; i++) {
        insn = &z->insns[i];
        switch (insn->op) {
        case GEMCA_RT_GUARD_BODY:
            guard_idx = insn->operand;
            break;
        case GEMCA_RT_PUSH_BODY:
        case GEMCA_RT_PUSH_VOXEL_BODY:
            if (sp >= OSH_GEMCA_RT_MAX_STACK) {
                return 0;
            }
            if (insn->operand >= 0 && (size_t) insn->operand < rt->nbodies) {
                stack[sp] = accel->body_aabb[insn->operand];
            } else {
                aabb_set_unbounded(&stack[sp]);
            }
            sp++;
            break;
        case GEMCA_RT_UNION:
            if (sp < 2) {
                return 0;
            }
            aabb_union(&stack[sp - 2], &stack[sp - 1]);
            sp--;
            break;
        case GEMCA_RT_INTERSECT:
            if (sp < 2) {
                return 0;
            }
            aabb_intersect(&stack[sp - 2], &stack[sp - 1]);
            sp--;
            break;
        case GEMCA_RT_DIFF:
            if (sp < 2) {
                return 0;
            }
            sp--; /* keep left operand */
            break;
        default:
            return 0;
        }
    }

    if (sp < 1) {
        return 0;
    }
    *out = stack[0];
    if (guard_idx >= 0 && (size_t) guard_idx < rt->nbodies) {
        aabb_intersect(out, &accel->body_aabb[guard_idx]);
    }
    return 1;
}

/* ---- Per-cell program specialization --------------------------------------- */

/* Stack entry for the constant folder: either a known constant 0 ("point is
 * never inside this subexpression anywhere in the cell") or an emitted
 * expression occupying out[start..] in the output buffer. */
struct fold_frame {
    int is_const0;
    int start;
};

/*
 * Specialize one zone program for one grid cell.  Bodies whose AABB does not
 * overlap the cell are folded to constant "outside"; boolean operators are
 * simplified, deleting instruction ranges that can no longer influence the
 * result.  Returns the specialized length, or -1 when the program folds to
 * constant false (the zone cannot contain any point of the cell).
 *
 * `out` must have room for z->ninsns entries.
 */
static int specialize_program(struct osh_gemca_accel const *accel,
                              struct gemca_rt_zone const *z,
                              struct gemca_rt_aabb const *cell_bb,
                              struct gemca_rt_insn *out) {
    struct fold_frame stack[OSH_GEMCA_RT_MAX_STACK];
    struct fold_frame a;
    struct fold_frame b;
    int sp = 0;
    int nout = 0;
    int i;
    int culled;
    struct gemca_rt_insn const *insn;

    for (i = 0; i < z->ninsns; i++) {
        insn = &z->insns[i];

        switch (insn->op) {
        case GEMCA_RT_GUARD_BODY:
            if (insn->operand >= 0 && accel->body_bounded[insn->operand]
                && !aabb_overlaps(&accel->body_aabb[insn->operand], cell_bb)) {
                /* Guard body never contains a cell point; the guard is an
                 * intersection term, so the whole zone is false here. */
                return -1;
            }
            out[nout].op = GEMCA_RT_GUARD_BODY;
            out[nout].operand = insn->operand;
            nout++;
            break;

        case GEMCA_RT_PUSH_BODY:
        case GEMCA_RT_PUSH_VOXEL_BODY:
            if (sp >= OSH_GEMCA_RT_MAX_STACK) {
                return -1;
            }
            culled =
                (insn->operand < 0)
                || (accel->body_bounded[insn->operand] && !aabb_overlaps(&accel->body_aabb[insn->operand], cell_bb));
            if (culled) {
                stack[sp].is_const0 = 1;
                stack[sp].start = nout;
            } else {
                stack[sp].is_const0 = 0;
                stack[sp].start = nout;
                out[nout] = *insn;
                nout++;
            }
            sp++;
            break;

        case GEMCA_RT_UNION:
        case GEMCA_RT_INTERSECT:
        case GEMCA_RT_DIFF:
            if (sp < 2) {
                return -1;
            }
            b = stack[sp - 1];
            a = stack[sp - 2];
            sp -= 2;

            if (insn->op == GEMCA_RT_UNION) {
                if (a.is_const0 && b.is_const0) {
                    stack[sp].is_const0 = 1;
                    stack[sp].start = a.start;
                } else if (a.is_const0) {
                    stack[sp] = b; /* 0 | x = x */
                } else if (b.is_const0) {
                    stack[sp] = a; /* x | 0 = x */
                } else {
                    out[nout].op = GEMCA_RT_UNION;
                    out[nout].operand = -1;
                    nout++;
                    stack[sp].is_const0 = 0;
                    stack[sp].start = a.start;
                }
            } else if (insn->op == GEMCA_RT_INTERSECT) {
                if (a.is_const0 || b.is_const0) {
                    nout = a.start; /* erase both operand programs */
                    stack[sp].is_const0 = 1;
                    stack[sp].start = a.start;
                } else {
                    out[nout].op = GEMCA_RT_INTERSECT;
                    out[nout].operand = -1;
                    nout++;
                    stack[sp].is_const0 = 0;
                    stack[sp].start = a.start;
                }
            } else { /* DIFF: a && !b */
                if (a.is_const0) {
                    nout = a.start;
                    stack[sp].is_const0 = 1;
                    stack[sp].start = a.start;
                } else if (b.is_const0) {
                    stack[sp] = a; /* x - 0 = x */
                } else {
                    out[nout].op = GEMCA_RT_DIFF;
                    out[nout].operand = -1;
                    nout++;
                    stack[sp].is_const0 = 0;
                    stack[sp].start = a.start;
                }
            }
            sp++;
            break;

        default:
            return -1;
        }
    }

    if (sp != 1 || stack[0].is_const0) {
        return -1;
    }
    return nout;
}

/* ---- Growable arrays -------------------------------------------------------- */

static int grow_u32(uint32_t **buf, size_t *cap, size_t need) {
    uint32_t *p;
    size_t newcap;
    if (need <= *cap) {
        return 1;
    }
    newcap = (*cap == 0) ? 256 : *cap;
    while (newcap < need) {
        newcap *= 2;
    }
    p = (uint32_t *) realloc(*buf, newcap * sizeof(uint32_t));
    if (!p) {
        return 0;
    }
    *buf = p;
    *cap = newcap;
    return 1;
}

static int grow_cands(struct gemca_accel_cand **buf, size_t *cap, size_t need) {
    struct gemca_accel_cand *p;
    size_t newcap;
    if (need <= *cap) {
        return 1;
    }
    newcap = (*cap == 0) ? 256 : *cap;
    while (newcap < need) {
        newcap *= 2;
    }
    p = (struct gemca_accel_cand *) realloc(*buf, newcap * sizeof(struct gemca_accel_cand));
    if (!p) {
        return 0;
    }
    *buf = p;
    *cap = newcap;
    return 1;
}

static int grow_insns(struct gemca_rt_insn **buf, size_t *cap, size_t need) {
    struct gemca_rt_insn *p;
    size_t newcap;
    if (need <= *cap) {
        return 1;
    }
    if (need > (size_t) ACCEL_MAX_POOL_INSNS) {
        return 0;
    }
    newcap = (*cap == 0) ? 1024 : *cap;
    while (newcap < need) {
        newcap *= 2;
    }
    p = (struct gemca_rt_insn *) realloc(*buf, newcap * sizeof(struct gemca_rt_insn));
    if (!p) {
        return 0;
    }
    *buf = p;
    *cap = newcap;
    return 1;
}

/* ---- Builder ---------------------------------------------------------------- */

void osh_gemca_accel_free(struct osh_gemca_accel *accel) {
    if (!accel) {
        return;
    }
    free(accel->cand_begin);
    free(accel->cands);
    free(accel->insns);
    free(accel->always);
    free(accel->body_aabb);
    free(accel->body_bounded);
    free(accel->zone_bounded);
    free(accel->zone_leaf_begin);
    free(accel->zone_leaves);
    free(accel->zone_leaf_unbounded);
    free(accel);
}

/* Compute the per-zone unique leaf-body lists used by the flat distance
 * evaluation.  Voxel zones are excluded (they dispatch to voxel traversal
 * before any RPN distance evaluation). */
static int build_zone_leaves(struct osh_gemca_runtime const *rt, struct osh_gemca_accel *accel) {
    uint8_t *seen;
    size_t iz;
    size_t total;
    int i;
    struct gemca_rt_zone const *z;

    accel->zone_leaf_begin = (uint32_t *) calloc(rt->nzones + 1, sizeof(uint32_t));
    if (!accel->zone_leaf_begin) {
        return 0;
    }
    seen = (uint8_t *) calloc(rt->nbodies, sizeof(uint8_t));
    if (!seen) {
        return 0;
    }

    total = 0;
    for (iz = 0; iz < rt->nzones; iz++) {
        z = &rt->zones[iz];
        accel->zone_leaf_begin[iz] = (uint32_t) total;
        if (z->voxel_body_idx >= 0) {
            continue;
        }
        for (i = 0; i < z->ninsns; i++) {
            if (z->insns[i].op == GEMCA_RT_PUSH_BODY && z->insns[i].operand >= 0
                && (size_t) z->insns[i].operand < rt->nbodies && !seen[z->insns[i].operand]) {
                seen[z->insns[i].operand] = 1;
                total++;
            }
        }
        for (i = 0; i < z->ninsns; i++) {
            if (z->insns[i].op == GEMCA_RT_PUSH_BODY && z->insns[i].operand >= 0
                && (size_t) z->insns[i].operand < rt->nbodies) {
                seen[z->insns[i].operand] = 0;
            }
        }
    }
    accel->zone_leaf_begin[rt->nzones] = (uint32_t) total;

    accel->zone_leaves = (uint32_t *) malloc((total ? total : 1) * sizeof(uint32_t));
    if (!accel->zone_leaves) {
        free(seen);
        return 0;
    }

    accel->zone_leaf_unbounded = (uint32_t *) calloc(rt->nzones, sizeof(uint32_t));
    if (!accel->zone_leaf_unbounded) {
        free(seen);
        return 0;
    }

    total = 0;
    for (iz = 0; iz < rt->nzones; iz++) {
        int pass;
        size_t zone_start;

        z = &rt->zones[iz];
        accel->zone_leaf_begin[iz] = (uint32_t) total;
        zone_start = total;
        if (z->voxel_body_idx >= 0) {
            continue;
        }
        /* Pass 0 collects unbounded leaves, pass 1 the bounded ones. */
        for (pass = 0; pass < 2; pass++) {
            for (i = 0; i < z->ninsns; i++) {
                if (z->insns[i].op == GEMCA_RT_PUSH_BODY && z->insns[i].operand >= 0
                    && (size_t) z->insns[i].operand < rt->nbodies && !seen[z->insns[i].operand]
                    && (int) accel->body_bounded[z->insns[i].operand] == pass) {
                    seen[z->insns[i].operand] = 1;
                    accel->zone_leaves[total++] = (uint32_t) z->insns[i].operand;
                }
            }
            if (pass == 0) {
                accel->zone_leaf_unbounded[iz] = (uint32_t) (total - zone_start);
            }
        }
        for (i = 0; i < z->ninsns; i++) {
            if (z->insns[i].op == GEMCA_RT_PUSH_BODY && z->insns[i].operand >= 0
                && (size_t) z->insns[i].operand < rt->nbodies) {
                seen[z->insns[i].operand] = 0;
            }
        }
    }
    accel->zone_leaf_begin[rt->nzones] = (uint32_t) total;

    free(seen);
    return 1;
}

enum osh_status osh_gemca_accel_build(struct osh_gemca_runtime *rt) {
    struct osh_gemca_accel *accel;
    struct gemca_rt_aabb *zone_bb = NULL;
    struct gemca_rt_insn *scratch = NULL;
    struct gemca_rt_aabb domain;
    struct gemca_rt_aabb cell_bb;
    double ext[3];
    double vol;
    double h;
    size_t target_cells;
    size_t iz;
    size_t ib;
    size_t ncells;
    size_t cell;
    size_t nbounded;
    size_t always_cap = 0;
    size_t cands_cap = 0;
    size_t insns_cap = 0;
    size_t ncands = 0;
    size_t ninsns_pool = 0;
    size_t max_zone_insns;
    uint32_t *cell_zone_count = NULL;
    uint32_t *cell_zone_begin = NULL;
    uint32_t *cell_zone_list = NULL;
    char const *env;
    int axis;
    int ix;
    int iy;
    int izc;
    int lo_idx[3];
    int hi_idx[3];
    int speclen;

    if (!rt) {
        return OSH_EINVAL;
    }
    rt->accel = NULL;

    env = getenv("OSH_GEMCA_ACCEL");
    if (env && env[0] == '0' && env[1] == '\0') {
        return OSH_OK;
    }
    if (rt->nzones == 0 || rt->nbodies == 0) {
        return OSH_OK;
    }

    /* Small geometries: the linear scan and plain RPN walks are already
     * cheap; indexing them costs more per query than it saves. */
    {
        size_t iz;
        size_t biggest = 0;
        for (iz = 0; iz < rt->nzones; iz++) {
            if ((size_t) rt->zones[iz].ninsns > biggest) {
                biggest = (size_t) rt->zones[iz].ninsns;
            }
        }
        if (rt->nzones < (size_t) OSH_GEMCA_ACCEL_MIN_ZONES && biggest < (size_t) OSH_GEMCA_ACCEL_MIN_ZONE_INSNS) {
            return OSH_OK;
        }
    }

    accel = (struct osh_gemca_accel *) calloc(1, sizeof(struct osh_gemca_accel));
    if (!accel) {
        return OSH_OK;
    }

    /* 1. Per-body universe AABBs. */
    accel->body_aabb = (struct gemca_rt_aabb *) calloc(rt->nbodies, sizeof(struct gemca_rt_aabb));
    accel->body_bounded = (uint8_t *) calloc(rt->nbodies, sizeof(uint8_t));
    if (!accel->body_aabb || !accel->body_bounded) {
        goto fail;
    }
    for (ib = 0; ib < rt->nbodies; ib++) {
        body_universe_aabb(rt, ib, &accel->body_aabb[ib], &accel->body_bounded[ib]);
    }

    /* 2. Per-zone unique leaf lists (for the flat distance evaluation). */
    if (!build_zone_leaves(rt, accel)) {
        goto fail;
    }

    /* 3. Per-zone AABBs and the grid domain. */
    zone_bb = (struct gemca_rt_aabb *) calloc(rt->nzones, sizeof(struct gemca_rt_aabb));
    accel->zone_bounded = (uint8_t *) calloc(rt->nzones, sizeof(uint8_t));
    if (!zone_bb || !accel->zone_bounded) {
        goto fail;
    }
    aabb_set_unbounded(&domain);
    domain.lo[0] = HUGE_VAL;
    domain.lo[1] = HUGE_VAL;
    domain.lo[2] = HUGE_VAL;
    domain.hi[0] = -HUGE_VAL;
    domain.hi[1] = -HUGE_VAL;
    domain.hi[2] = -HUGE_VAL;
    nbounded = 0;
    max_zone_insns = 0;
    for (iz = 0; iz < rt->nzones; iz++) {
        if ((size_t) rt->zones[iz].ninsns > max_zone_insns) {
            max_zone_insns = (size_t) rt->zones[iz].ninsns;
        }
        if (zone_aabb(rt, accel, &rt->zones[iz], &zone_bb[iz]) && aabb_is_bounded(&zone_bb[iz])) {
            accel->zone_bounded[iz] = 1;
            aabb_union(&domain, &zone_bb[iz]);
            nbounded++;
        } else {
            accel->zone_bounded[iz] = 0;
            if (!grow_u32(&accel->always, &always_cap, accel->nalways + 1)) {
                goto fail;
            }
            accel->always[accel->nalways++] = (uint32_t) iz;
        }
    }
    if (nbounded == 0 || !aabb_is_bounded(&domain)) {
        goto fail; /* nothing to index; fall back to linear scan */
    }

    /* 4. Grid resolution: roughly cubic cells, ACCEL_CELLS_PER_ZONE cells
     *    per zone, clamped per axis. */
    for (axis = 0; axis < 3; axis++) {
        ext[axis] = domain.hi[axis] - domain.lo[axis];
        if (!(ext[axis] > 0.0)) {
            ext[axis] = 1e-6;
        }
    }
    target_cells = rt->nzones * (size_t) ACCEL_CELLS_PER_ZONE;
    if (target_cells < (size_t) ACCEL_MIN_CELLS) {
        target_cells = (size_t) ACCEL_MIN_CELLS;
    }
    vol = ext[0] * ext[1] * ext[2];
    h = cbrt(vol / (double) target_cells);
    ncells = 1;
    for (axis = 0; axis < 3; axis++) {
        int na = (int) ceil(ext[axis] / h);
        if (na < 1) {
            na = 1;
        }
        if (na > ACCEL_MAX_CELLS_PER_AXIS) {
            na = ACCEL_MAX_CELLS_PER_AXIS;
        }
        accel->n[axis] = na;
        accel->lo[axis] = domain.lo[axis];
        accel->hi[axis] = domain.hi[axis];
        accel->inv_cell[axis] = (double) na / ext[axis];
        ncells *= (size_t) na;
    }
    accel->ncells = ncells;

    /* 5. Rasterize bounded zone AABBs into per-cell zone lists (CSR).
     *    Zones are visited in ascending index order, so each cell list is
     *    sorted by construction — required to preserve first-match-wins. */
    cell_zone_count = (uint32_t *) calloc(ncells + 1, sizeof(uint32_t));
    if (!cell_zone_count) {
        goto fail;
    }
    for (iz = 0; iz < rt->nzones; iz++) {
        if (!accel->zone_bounded[iz]) {
            continue;
        }
        for (axis = 0; axis < 3; axis++) {
            lo_idx[axis] = (int) floor((zone_bb[iz].lo[axis] - accel->lo[axis]) * accel->inv_cell[axis]);
            hi_idx[axis] = (int) floor((zone_bb[iz].hi[axis] - accel->lo[axis]) * accel->inv_cell[axis]);
            if (lo_idx[axis] < 0) {
                lo_idx[axis] = 0;
            }
            if (hi_idx[axis] >= accel->n[axis]) {
                hi_idx[axis] = accel->n[axis] - 1;
            }
        }
        for (izc = lo_idx[2]; izc <= hi_idx[2]; izc++) {
            for (iy = lo_idx[1]; iy <= hi_idx[1]; iy++) {
                for (ix = lo_idx[0]; ix <= hi_idx[0]; ix++) {
                    cell = ((size_t) izc * (size_t) accel->n[1] + (size_t) iy) * (size_t) accel->n[0] + (size_t) ix;
                    cell_zone_count[cell]++;
                }
            }
        }
    }
    cell_zone_begin = (uint32_t *) malloc((ncells + 1) * sizeof(uint32_t));
    if (!cell_zone_begin) {
        goto fail;
    }
    cell_zone_begin[0] = 0;
    for (cell = 0; cell < ncells; cell++) {
        cell_zone_begin[cell + 1] = cell_zone_begin[cell] + cell_zone_count[cell];
        cell_zone_count[cell] = 0;
    }
    cell_zone_list = (uint32_t *) malloc((cell_zone_begin[ncells] ? cell_zone_begin[ncells] : 1) * sizeof(uint32_t));
    if (!cell_zone_list) {
        goto fail;
    }
    for (iz = 0; iz < rt->nzones; iz++) {
        if (!accel->zone_bounded[iz]) {
            continue;
        }
        for (axis = 0; axis < 3; axis++) {
            lo_idx[axis] = (int) floor((zone_bb[iz].lo[axis] - accel->lo[axis]) * accel->inv_cell[axis]);
            hi_idx[axis] = (int) floor((zone_bb[iz].hi[axis] - accel->lo[axis]) * accel->inv_cell[axis]);
            if (lo_idx[axis] < 0) {
                lo_idx[axis] = 0;
            }
            if (hi_idx[axis] >= accel->n[axis]) {
                hi_idx[axis] = accel->n[axis] - 1;
            }
        }
        for (izc = lo_idx[2]; izc <= hi_idx[2]; izc++) {
            for (iy = lo_idx[1]; iy <= hi_idx[1]; iy++) {
                for (ix = lo_idx[0]; ix <= hi_idx[0]; ix++) {
                    cell = ((size_t) izc * (size_t) accel->n[1] + (size_t) iy) * (size_t) accel->n[0] + (size_t) ix;
                    cell_zone_list[cell_zone_begin[cell] + cell_zone_count[cell]] = (uint32_t) iz;
                    cell_zone_count[cell]++;
                }
            }
        }
    }

    /* 6. Specialize each (cell, zone) program by constant folding. */
    scratch = (struct gemca_rt_insn *) malloc((max_zone_insns ? max_zone_insns : 1) * sizeof(struct gemca_rt_insn));
    accel->cand_begin = (uint32_t *) malloc((ncells + 1) * sizeof(uint32_t));
    if (!scratch || !accel->cand_begin) {
        goto fail;
    }
    for (cell = 0; cell < ncells; cell++) {
        size_t k;

        accel->cand_begin[cell] = (uint32_t) ncands;

        ix = (int) (cell % (size_t) accel->n[0]);
        iy = (int) ((cell / (size_t) accel->n[0]) % (size_t) accel->n[1]);
        izc = (int) (cell / ((size_t) accel->n[0] * (size_t) accel->n[1]));
        cell_bb.lo[0] = accel->lo[0] + (double) ix / accel->inv_cell[0];
        cell_bb.hi[0] = accel->lo[0] + (double) (ix + 1) / accel->inv_cell[0];
        cell_bb.lo[1] = accel->lo[1] + (double) iy / accel->inv_cell[1];
        cell_bb.hi[1] = accel->lo[1] + (double) (iy + 1) / accel->inv_cell[1];
        cell_bb.lo[2] = accel->lo[2] + (double) izc / accel->inv_cell[2];
        cell_bb.hi[2] = accel->lo[2] + (double) (izc + 1) / accel->inv_cell[2];

        for (k = cell_zone_begin[cell]; k < cell_zone_begin[cell + 1]; k++) {
            uint32_t zidx = cell_zone_list[k];

            speclen = specialize_program(accel, &rt->zones[zidx], &cell_bb, scratch);
            if (speclen < 0) {
                continue; /* constant false: zone unreachable from this cell */
            }
            if (!grow_insns(&accel->insns, &insns_cap, ninsns_pool + (size_t) speclen)
                || !grow_cands(&accel->cands, &cands_cap, ncands + 1)) {
                goto fail;
            }
            memcpy(&accel->insns[ninsns_pool], scratch, (size_t) speclen * sizeof(struct gemca_rt_insn));
            accel->cands[ncands].zone_idx = zidx;
            accel->cands[ncands].insn_begin = (uint32_t) ninsns_pool;
            accel->cands[ncands].ninsns = (uint32_t) speclen;
            ninsns_pool += (size_t) speclen;
            ncands++;
        }
    }
    accel->cand_begin[ncells] = (uint32_t) ncands;

    free(scratch);
    free(cell_zone_count);
    free(cell_zone_begin);
    free(cell_zone_list);
    free(zone_bb);

    rt->accel = accel;
    return OSH_OK;

fail:
    free(scratch);
    free(cell_zone_count);
    free(cell_zone_begin);
    free(cell_zone_list);
    free(zone_bb);
    osh_gemca_accel_free(accel);
    return OSH_OK;
}
