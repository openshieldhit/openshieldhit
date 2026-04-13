#include "gemca/runtime/osh_gemca_runtime.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_coord.h"
#include "common/osh_logger.h"
#include "common/osh_ray.h"
#include "common/osh_vect.h"
#include "gemca/osh_gemca2.h"
#include "gemca/osh_gemca2_defines.h"

/* ---- Local constants ----------------------------------------------------- */

/*
 * Convenience macro for selecting the smaller positive value of two doubles.
 * Identical in purpose to the MIN macro in osh_gemca2_dist.c but scoped to
 * this translation unit.
 */
#define _RT_MIN(a, b) (((a) < (b)) ? (a) : (b))

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
static enum osh_status setup_surfaces(struct gemca_workspace const *wg, struct gemca_runtime *rt);
static enum osh_status setup_bodies(struct gemca_workspace const *wg, struct gemca_runtime *rt);
static enum osh_status setup_zones(struct gemca_workspace const *wg, struct gemca_runtime *rt);

/* Zone compilation */
static enum osh_status compile_zone(struct zone const *z,
                                    struct gemca_workspace const *wg,
                                    struct gemca_rt_zone *zrt);
static void compile_node(struct cgnode const *node,
                         struct gemca_workspace const *wg,
                         struct gemca_rt_insn *insns,
                         int *ninsns);
static int count_leaves(struct cgnode const *node);
static struct body const *find_guard_body(struct cgnode const *node);
static int find_body_index(struct gemca_workspace const *wg, struct body const *b);

/* RPN evaluators */
static int eval_membership(struct gemca_runtime const *rt,
                           struct gemca_rt_zone const *z,
                           struct ray const *r);
static double eval_distance(struct gemca_runtime const *rt,
                            struct gemca_rt_zone const *z,
                            struct ray const *r,
                            int *is_inside);

/* Body evaluators */
static int in_body_rt(struct gemca_runtime const *rt, int body_idx, struct ray const *r);
static double dist_body_rt(struct gemca_runtime const *rt, int body_idx, struct ray const *r);

/* Ray transform */
static enum osh_status transform_to_local_rt(struct gemca_rt_body const *b,
                                             struct ray const *r,
                                             struct ray *tr);

/* Surface evaluators */
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
enum osh_status osh_gemca_runtime_setup(struct gemca_workspace const *wg, struct gemca_runtime *rt) {
    enum osh_status rc;

    if (!wg || !rt) {
        return OSH_EINVAL;
    }

    rt->workspace = wg;

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

    rc = setup_zones(wg, rt);
    if (rc != OSH_OK) {
        osh_gemca_runtime_free(rt);
        return rc;
    }

    return OSH_OK;
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
void osh_gemca_runtime_free(struct gemca_runtime *rt) {
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
size_t osh_gemca_runtime_get_zone(struct gemca_runtime const *rt, struct ray const *r) {
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
double osh_gemca_runtime_get_distance(struct gemca_runtime const *rt, size_t zone_idx, struct ray const *r) {
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
        if (d < 0.0) {
            osh_error("osh_gemca_runtime_get_distance(): negative step distance");
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
static enum osh_status setup_surfaces(struct gemca_workspace const *wg, struct gemca_runtime *rt) {
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
 * Records the transformation matrix, coordinate system, body type, and the
 * offset (surf_begin) into rt->surfaces[] for each body.  The offset is
 * computed by accumulating surface counts in body order, which matches the
 * layout produced by setup_surfaces().
 *
 * @param[in]  wg  Cold workspace.
 * @param[out] rt  Runtime (rt->surfaces must already be populated).
 *
 * @returns OSH_OK or OSH_ENOMEM.
 */
static enum osh_status setup_bodies(struct gemca_workspace const *wg, struct gemca_runtime *rt) {
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
static enum osh_status setup_zones(struct gemca_workspace const *wg, struct gemca_runtime *rt) {
    size_t iz;
    enum osh_status rc;

    rt->zones = (struct gemca_rt_zone *) calloc(wg->nzones, sizeof(struct gemca_rt_zone));
    if (!rt->zones) {
        return OSH_ENOMEM;
    }
    rt->nzones = wg->nzones;

    for (iz = 0; iz < wg->nzones; iz++) {
        rc = compile_zone(wg->zones[iz], wg, &rt->zones[iz]);
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
                                    struct gemca_workspace const *wg,
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

    /* Prepend guard if a suitable body can be identified. */
    guard = find_guard_body(&z->node);
    if (guard) {
        guard_idx = find_body_index(wg, guard);
        if (guard_idx >= 0) {
            zrt->insns[ninsns].op = GEMCA_RT_GUARD_BODY;
            zrt->insns[ninsns].operand = guard_idx;
            ninsns++;
        }
    }

    /* Compile the CSG tree post-order into the remaining slots. */
    compile_node(&z->node, wg, zrt->insns, &ninsns);

    zrt->ninsns = ninsns;
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
                         struct gemca_workspace const *wg,
                         struct gemca_rt_insn *insns,
                         int *ninsns) {
    int body_idx;

    if (node->type == _OSH_GEMCA_CGNODE_BODY) {
        body_idx = find_body_index(wg, node->body);
        insns[*ninsns].op = (node->body->type == OSH_GEMCA_BODY_VOX)
                                ? GEMCA_RT_PUSH_VOXEL_BODY
                                : GEMCA_RT_PUSH_BODY;
        insns[*ninsns].operand = body_idx;
        (*ninsns)++;
        return;
    }

    /* Interior node: left subtree, right subtree, then operator. */
    compile_node(node->left, wg, insns, ninsns);
    compile_node(node->right, wg, insns, ninsns);

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
        osh_error("compile_node(): unknown CSG operator '%c'", node->op);
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
static int find_body_index(struct gemca_workspace const *wg, struct body const *b) {
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
static int eval_membership(struct gemca_runtime const *rt,
                           struct gemca_rt_zone const *z,
                           struct ray const *r) {
    int stack[OSH_GEMCA_RT_MAX_STACK];
    int sp;
    int i;
    struct gemca_rt_insn const *insn;

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
            stack[sp++] = in_body_rt(rt, insn->operand, r);
            break;

        case GEMCA_RT_UNION:
            stack[sp - 2] = stack[sp - 2] || stack[sp - 1];
            sp--;
            break;

        case GEMCA_RT_INTERSECT:
            stack[sp - 2] = stack[sp - 2] && stack[sp - 1];
            sp--;
            break;

        case GEMCA_RT_DIFF:
            /* left && !right  (left = sp-2, right = sp-1 in RPN order) */
            stack[sp - 2] = stack[sp - 2] && !stack[sp - 1];
            sp--;
            break;

        default:
            osh_error("eval_membership(): unknown opcode %d", insn->op);
            return 0;
        }
    }

    return (sp > 0) ? stack[0] : 0;
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
 */
static double eval_distance(struct gemca_runtime const *rt,
                            struct gemca_rt_zone const *z,
                            struct ray const *r,
                            int *is_inside) {
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
        case GEMCA_RT_PUSH_VOXEL_BODY:
            /*
             * TODO: for PUSH_VOXEL_BODY and z->voxel_body_idx >= 0, dispatch to
             * the Jacobs voxel traversal algorithm here instead of the RPP
             * boundary distance.  Current fallback: treat as a regular body.
             */
            stack[sp].is_inside = in_body_rt(rt, insn->operand, r);
            stack[sp].dist = dist_body_rt(rt, insn->operand, r);
            sp++;
            break;

        case GEMCA_RT_UNION:
            b = stack[--sp];
            a = stack[sp - 1];
            stack[sp - 1].is_inside = a.is_inside || b.is_inside;
            stack[sp - 1].dist = _minpos(a.dist, b.dist);
            break;

        case GEMCA_RT_INTERSECT:
            b = stack[--sp];
            a = stack[sp - 1];
            stack[sp - 1].is_inside = a.is_inside && b.is_inside;
            stack[sp - 1].dist = _minpos(a.dist, b.dist);
            break;

        case GEMCA_RT_DIFF:
            b = stack[--sp];
            a = stack[sp - 1];
            stack[sp - 1].is_inside = a.is_inside && !b.is_inside;
            stack[sp - 1].dist = _minpos(a.dist, b.dist);
            break;

        default:
            osh_error("eval_distance(): unknown opcode %d", insn->op);
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
static int in_body_rt(struct gemca_runtime const *rt, int body_idx, struct ray const *r) {
    struct gemca_rt_body const *b;
    struct ray tr;
    int i;

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
static double dist_body_rt(struct gemca_runtime const *rt, int body_idx, struct ray const *r) {
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
static enum osh_status transform_to_local_rt(struct gemca_rt_body const *b,
                                             struct ray const *r,
                                             struct ray *tr) {
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
        osh_error("transform_to_local_rt(): unsupported coordinate system %d", (int) b->coord);
        return OSH_ENOTSUP;
    }

    return OSH_OK;
}

/* ---- Surface evaluators -------------------------------------------------- */

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
    double d;
    double dot;
    int i;

    switch (sf->type) {

    case OSH_GEMCA_SURF_SPHERE:
        d = osh_vect_len2(r->p) - sf->p[0];
        if (d > OSH_GEMCA_SMALL) {
            return 0;
        } else if (d < -OSH_GEMCA_SMALL) {
            return 1;
        }
        return (osh_vect_dot(r->p, r->cp) < 0.0) ? 1 : 0;

    case OSH_GEMCA_SURF_ELLIPSOID:
        d = 0.0;
        for (i = 0; i < 3; i++) {
            d += (r->p[i] * r->p[i]) / sf->p[i];
        }
        if (d > 1.0 + OSH_GEMCA_SMALL) {
            return 0;
        } else if (d < 1.0 - OSH_GEMCA_SMALL) {
            return 1;
        }
        dot = 0.0;
        for (i = 0; i < 3; i++) {
            dot += (r->p[i] / sf->p[i]) * r->cp[i];
        }
        return (dot < 0.0) ? 1 : 0;

    case OSH_GEMCA_SURF_CYLZ:
        d = r->p[0] * r->p[0] + r->p[1] * r->p[1] - sf->p[0];
        if (d > OSH_GEMCA_SMALL) {
            return 0;
        } else if (d < -OSH_GEMCA_SMALL) {
            return 1;
        }
        dot = r->p[0] * r->cp[0] + r->p[1] * r->cp[1];
        return (dot < 0.0) ? 1 : 0;

    case OSH_GEMCA_SURF_ELLZ:
        d = (r->p[0] * r->p[0] / sf->p[0]) + (r->p[1] * r->p[1] / sf->p[1]) - 1.0;
        if (d > OSH_GEMCA_SMALL) {
            return 0;
        } else if (d < -OSH_GEMCA_SMALL) {
            return 1;
        }
        dot = (r->p[0] / sf->p[0]) * r->cp[0] + (r->p[1] / sf->p[1]) * r->cp[1];
        return (dot < 0.0) ? 1 : 0;

    case OSH_GEMCA_SURF_CONE:
        d = (r->p[0] * r->p[0]) + (r->p[1] * r->p[1]) - sf->p[1] * (r->p[2] * r->p[2]);
        if (d > OSH_GEMCA_SMALL) {
            return 0;
        } else if (d < -OSH_GEMCA_SMALL) {
            return 1;
        }
        dot = (r->p[0] * r->cp[0]) + (r->p[1] * r->cp[1]) - (sf->p[1] * r->p[2] * r->cp[2]);
        return (dot < 0.0) ? 1 : 0;

    case OSH_GEMCA_SURF_PLANEX:
        d = sf->p[0] * r->p[0] + sf->p[1];
        if (d > OSH_GEMCA_SMALL) {
            return 0;
        } else if (d < -OSH_GEMCA_SMALL) {
            return 1;
        }
        return (sf->p[0] * r->cp[0] > 0.0) ? 0 : 1;

    case OSH_GEMCA_SURF_PLANEY:
        d = sf->p[0] * r->p[1] + sf->p[1];
        if (d > OSH_GEMCA_SMALL) {
            return 0;
        } else if (d < -OSH_GEMCA_SMALL) {
            return 1;
        }
        return (sf->p[0] * r->cp[1] > 0.0) ? 0 : 1;

    case OSH_GEMCA_SURF_PLANEZ:
        d = sf->p[0] * r->p[2] + sf->p[1];
        if (d > OSH_GEMCA_SMALL) {
            return 0;
        } else if (d < -OSH_GEMCA_SMALL) {
            return 1;
        }
        return (sf->p[0] * r->cp[2] > 0.0) ? 0 : 1;

    case OSH_GEMCA_SURF_PLANE:
        /* p[0..2] = normal (A,B,C); p[3] = offset D.  Ax+By+Cz+D <= 0 is inside. */
        d = osh_vect_dot(sf->p, r->p) + sf->p[3];
        if (d > OSH_GEMCA_SMALL) {
            return 0;
        } else if (d < -OSH_GEMCA_SMALL) {
            return 1;
        }
        return (osh_vect_dot(sf->p, r->cp) > 0.0) ? 0 : 1;

    default:
        osh_error("_check_surface_rt(): unknown surface type %d", sf->type);
        return 1;
    }
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
 * non-positive.  The tangential case (discriminant == 0) is treated as no
 * intersection.
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
