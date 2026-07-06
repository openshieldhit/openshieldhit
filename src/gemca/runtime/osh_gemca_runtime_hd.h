#ifndef OSH_GEMCA_RUNTIME_HD_H
#define OSH_GEMCA_RUNTIME_HD_H

/*
 * osh_gemca_runtime_hd.h — device-compilable zone-membership evaluators.
 *
 * Function bodies are moved verbatim from osh_gemca_runtime.c and marked
 * OSH_HD static inline so the identical lines compile both as host code
 * (osh_gemca_runtime.c includes this header and delegates to these twins)
 * and as CUDA device code (kernels in src/gpu/ call them directly).
 *
 * There is deliberately NO independent GPU implementation of any of this
 * logic: the coordinate-system dispatch in the body transform (only the
 * translation column of t[] is populated for OSH_COORD_BCALIGN bodies),
 * the on-surface direction disambiguation, and the RPN operand order are
 * all easy to get subtly wrong when re-derived by hand.  Single source is
 * the correctness strategy, not an optimisation.
 *
 * Signature convention: the twins take the runtime's flat arrays
 * (surfaces/bodies + counts) instead of struct osh_gemca_runtime, because
 * device code receives those arrays as device pointers inside a POD view
 * struct and never sees the host runtime object.
 */

#include "common/osh_coord.h"
#include "common/osh_hd.h"
#include "common/osh_ray.h"
#include "common/osh_ray_hd.h"
#include "gemca/osh_gemca2.h"
#include "gemca/osh_gemca2_defines.h"
#include "gemca/runtime/osh_gemca_runtime.h"

/*
 * _osh_gemca_rt_transform_to_local_hd() — body of transform_to_local_rt()
 * (see osh_gemca_runtime.c for the full contract).  Transforms a ray from
 * OSH_COORD_UNIVERSE into a body's local coordinate system:
 *
 *   OSH_COORD_UNIVERSE — identity; ray is copied unchanged.
 *   OSH_COORD_BCALIGN  — pure translation via the last column of t[]
 *                        (the rotation block of t[] is NOT populated).
 *   OSH_COORD_BZALIGN  — full affine transform via the ray transform.
 */
OSH_HD static inline enum osh_status
_osh_gemca_rt_transform_to_local_hd(struct gemca_rt_body const *b, struct ray const *r, struct ray *tr) {
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
        _osh_ray_transform_hd(r, tr, b->t);
        break;

    default:
        return OSH_ENOTSUP;
    }

    return OSH_OK;
}

/*
 * _osh_gemca_rt_check_surface_components_hd() — body of
 * _check_surface_components_rt() (see osh_gemca_runtime.c).  Inside check
 * for one flat surface against a body-local position/direction; a negative
 * value of the surface implicit function means "inside"; on-surface cases
 * are disambiguated by the ray direction.
 */
OSH_HD static inline int _osh_gemca_rt_check_surface_components_hd(
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

/*
 * _osh_gemca_rt_check_surface_hd() — body of _check_surface_rt().
 * Ray already transformed to body-local coordinates.
 */
OSH_HD static inline int _osh_gemca_rt_check_surface_hd(struct gemca_rt_surface const *sf, struct ray const *r) {
    return _osh_gemca_rt_check_surface_components_hd(sf, r->p[0], r->p[1], r->p[2], r->cp[0], r->cp[1], r->cp[2]);
}

/*
 * _osh_gemca_rt_in_body_hd() — body of in_body_rt().  A ray is inside a
 * body iff it is inside every surface of the body's flat surface slice.
 */
OSH_HD static inline int _osh_gemca_rt_in_body_hd(struct gemca_rt_surface const *surfaces,
                                                  struct gemca_rt_body const *bodies,
                                                  size_t nbodies,
                                                  int body_idx,
                                                  struct ray const *r) {
    struct gemca_rt_body const *b;
    struct ray tr;
    int i;

    if (body_idx < 0 || (size_t) body_idx >= nbodies) {
        return 0;
    }

    b = &bodies[body_idx];

    if (_osh_gemca_rt_transform_to_local_hd(b, r, &tr) != OSH_OK) {
        return 0;
    }

    for (i = 0; i < b->nsurfs; i++) {
        if (!_osh_gemca_rt_check_surface_hd(&surfaces[b->surf_begin + (size_t) i], &tr)) {
            return 0;
        }
    }
    return 1;
}

/*
 * _osh_gemca_rt_eval_membership_hd() — body of eval_membership().
 * Evaluates one zone's RPN instruction array (either the per-zone insns[]
 * on the host or a slice of insns_flat[] on the device) left-to-right with
 * a small integer stack.  GUARD_BODY at the head fast-rejects; PUSH leaves
 * a membership flag; UNION/INTERSECT/DIFF fold the top two entries.
 */
OSH_HD static inline int _osh_gemca_rt_eval_membership_hd(struct gemca_rt_surface const *surfaces,
                                                          struct gemca_rt_body const *bodies,
                                                          size_t nbodies,
                                                          struct gemca_rt_insn const *insns,
                                                          int ninsns,
                                                          struct ray const *r) {
    int stack[OSH_GEMCA_RT_MAX_STACK];
    int sp;
    int i;
    struct gemca_rt_insn const *insn;

    /*
     * Fast path: single-instruction zones (one PUSH_BODY after the guard
     * optimisation removed guards for leaf zones) need no stack machinery.
     * This is the common case for simple body-per-zone geometry.
     */
    if (ninsns == 1) {
        return _osh_gemca_rt_in_body_hd(surfaces, bodies, nbodies, insns[0].operand, r);
    }

    sp = 0;

    for (i = 0; i < ninsns; i++) {
        insn = &insns[i];

        switch (insn->op) {

        case GEMCA_RT_GUARD_BODY:
            /* Fast-reject: if outside guard body, skip the entire zone. */
            if (!_osh_gemca_rt_in_body_hd(surfaces, bodies, nbodies, insn->operand, r)) {
                return 0;
            }
            break;

        case GEMCA_RT_PUSH_BODY:
        case GEMCA_RT_PUSH_VOXEL_BODY:
            if (sp >= OSH_GEMCA_RT_MAX_STACK) {
                return 0;
            }
            stack[sp++] = _osh_gemca_rt_in_body_hd(surfaces, bodies, nbodies, insn->operand, r);
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

/*
 * _osh_gemca_rt_find_zone_flat_hd() — zone lookup over the flat instruction
 * store (osh_gemca_runtime.insns_flat / insn_begin).  Semantics mirror
 * osh_gemca_runtime_get_zone(): first zone whose expression evaluates true,
 * or OSH_GEMCA_ZONE_INDEX_INVALID.  The flat store carries the same
 * instructions as the per-zone arrays (unit-tested equivalence), so this
 * walks the identical evaluator over identical data.
 */
OSH_HD static inline size_t _osh_gemca_rt_find_zone_flat_hd(struct gemca_rt_surface const *surfaces,
                                                            struct gemca_rt_body const *bodies,
                                                            size_t nbodies,
                                                            struct gemca_rt_insn const *insns_flat,
                                                            int const *insn_begin,
                                                            size_t nzones,
                                                            struct ray const *r) {
    size_t i;
    int begin;
    int ninsns;

    for (i = 0; i < nzones; i++) {
        begin = insn_begin[i];
        ninsns = insn_begin[i + 1] - begin;
        if (_osh_gemca_rt_eval_membership_hd(surfaces, bodies, nbodies, &insns_flat[begin], ninsns, r)) {
            return i;
        }
    }

    return OSH_GEMCA_ZONE_INDEX_INVALID;
}

#endif /* OSH_GEMCA_RUNTIME_HD_H */
