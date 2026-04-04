#include "beam/osh_beam_model.h"

#include <math.h>

#include "common/osh_const.h"
#include "common/osh_coord.h"
#include "common/osh_rc.h"
#include "common/osh_vect.h"

/* ---- Step 1: spot selection ---------------------------------------------- */

/* Return a pointer to the spot to use for this history.
 * Single-spot and phase-space modes always return spots[0].
 * SOBP mode draws a spot from the list weighted by spot->wt.
 * TODO: implement weighted draw for SOBP. */
static struct beam_spot const *_select_spot(struct beam_workspace const *wb) {
    /* TODO: SOBP weighted selection */
    return &wb->spots[0];
}

/* ---- Step 2: energy sampling --------------------------------------------- */

/* Sample total kinetic energy [MeV] from the spot's energy distribution.
 * Current implementation uses Gaussian sampling around t0 with sigma tsigma,
 * clamped at zero. */
static double _sample_energy(struct osh_rng *rng, struct beam_spot const *spot) {
    double e;

    if (spot->tsigma <= 0.0) {
        return spot->t0;
    }

    e = osh_rng_gauss(rng, spot->t0, spot->tsigma);
    if (e < 0.0) {
        e = 0.0;
    }
    return e;
}

/* ---- Step 3: transverse phase-space sampling ----------------------------- */

static double _sample_uniform_symmetric(struct osh_rng *rng, double half_width) {
    return (2.0 * osh_rng_double(rng) - 1.0) * half_width;
}

static double _shape_sigma_1d(struct beam_spot const *spot, int axis) {
    double rin;
    double rout;

    switch (spot->shape) {
    case OSH_BEAM_SHAPE_PENCIL:
        return 0.0;

    case OSH_BEAM_SHAPE_GAUSSIAN:
        return spot->size[axis];

    case OSH_BEAM_SHAPE_SQUARE:
        return spot->size[axis] / sqrt(3.0);

    case OSH_BEAM_SHAPE_CIRCULAR:
        rin = fmin(spot->size[0], spot->size[1]);
        rout = fmax(spot->size[0], spot->size[1]);
        return 0.5 * sqrt(rin * rin + rout * rout);

    default:
        return -1.0;
    }
}

static int _sample_position(struct osh_rng *rng, struct beam_spot const *spot, double *x_out, double *y_out) {
    double rin;
    double rout;
    double r2;
    double phi;

    switch (spot->shape) {
    case OSH_BEAM_SHAPE_PENCIL:
        *x_out = 0.0;
        *y_out = 0.0;
        return OSH_OK;

    case OSH_BEAM_SHAPE_GAUSSIAN:
        *x_out = (spot->size[0] > 0.0) ? spot->size[0] * osh_rng_gauss01(rng) : 0.0;
        *y_out = (spot->size[1] > 0.0) ? spot->size[1] * osh_rng_gauss01(rng) : 0.0;
        return OSH_OK;

    case OSH_BEAM_SHAPE_SQUARE:
        *x_out = _sample_uniform_symmetric(rng, spot->size[0]);
        *y_out = _sample_uniform_symmetric(rng, spot->size[1]);
        return OSH_OK;

    case OSH_BEAM_SHAPE_CIRCULAR:
        rin = fmin(spot->size[0], spot->size[1]);
        rout = fmax(spot->size[0], spot->size[1]);
        if (rin < 0.0 || rout < rin) {
            return OSH_EINVAL;
        }
        phi = 2.0 * OSH_M_PI * osh_rng_double(rng);
        r2 = rin * rin + (rout * rout - rin * rin) * osh_rng_double(rng);
        rout = sqrt(r2);
        *x_out = rout * cos(phi);
        *y_out = rout * sin(phi);
        return OSH_OK;

    default:
        return OSH_EINVAL;
    }
}

/* Sample Gaussian angle in one transverse plane using the Fermi-Eyges
 * correlation form. sigma_x and sigma_xp are the RMS widths; rho is the
 * correlation coefficient. u2 is an independent N(0,1) sample.
 *
 *   x  = sigma_x  * u1
 *   x' = sigma_xp * (rho * u1 + sqrt(1 - rho^2) * u2)
 *
 * For the current single-spot path we use the equivalent closed form
 * directly. For non-Gaussian shapes we feed in the actual sampled x and the
 * RMS width corresponding to that shape, which preserves the requested
 * position-angle correlation. */
static void _sample_angle_1d(
    struct osh_rng *rng, double x, double sigma_x, double sigma_xp, double rho, double *xp_out) {
    double u2;
    double rho_clamped;
    double rho_orth;

    if (sigma_xp <= 0.0) {
        *xp_out = 0.0;
        return;
    }

    if (sigma_x <= 0.0) {
        *xp_out = sigma_xp * osh_rng_gauss01(rng);
        return;
    }

    rho_clamped = rho;
    if (rho_clamped > 1.0) {
        rho_clamped = 1.0;
    } else if (rho_clamped < -1.0) {
        rho_clamped = -1.0;
    }
    rho_orth = sqrt(fmax(0.0, 1.0 - rho_clamped * rho_clamped));
    u2 = osh_rng_gauss01(rng);
    *xp_out = sigma_xp * (rho_clamped * (x / sigma_x) + rho_orth * u2);
}

/* ---- Step 4: SAD correction ---------------------------------------------- */

/* Apply source-axis-distance correction in PZALIGN.
 *
 * SAD is a positive source-to-isocenter distance to an upstream virtual
 * source, not a signed z coordinate. It is defined relative to
 * isocenter/focal geometry, not relative to the beam start plane. The actual
 * local start position is therefore
 *   spot->p + ray->p
 * and the fan-out slope for each plane is taken from a virtual source at
 *   z = -sad[axis]
 * pointing to that start point. Since BEAMPOS is already backprojected by the
 * caller, this updates only the direction, not the start position. */
static int _apply_sad(struct ray_v *ray, struct beam_spot const *spot, struct beam_shared const *sh) {
    double x;
    double y;
    double z;
    double dzx;
    double dzy;

    x = spot->p[0] + ray->p[0];
    y = spot->p[1] + ray->p[1];
    z = spot->p[2] + ray->p[2];

    dzx = z + sh->sad[0];
    dzy = z + sh->sad[1];
    if (dzx <= 0.0 || dzy <= 0.0) {
        return OSH_ESTATE;
    }

    ray->v[0] += x / dzx;
    ray->v[1] += y / dzy;
    return OSH_OK;
}

/* ---- Step 5: standard affine transform PZALIGN → UNIVERSE ---------------- */

/* Apply the beam-model sampling matrix:
 *   p_u = R * p_l + t
 *   v_u = R * v_l
 * where R,t are stored in spot->_tm. This intentionally does NOT use
 * osh_coord_trans_ray(), which follows the legacy SHIELD-HIT/GEMCA sign
 * convention for translation. */
static void _apply_transform(struct ray_v *ray, struct beam_spot const *spot) {
    double p[3];
    double v[3];

    osh_vect_copy(ray->p, p);
    osh_vect_copy(ray->v, v);

    osh_vect_trans_point_affine(p, ray->p, spot->_tm);
    osh_vect_trans_vector_affine(v, ray->v, spot->_tm);

    ray->system = OSH_COORD_UNIVERSE;
}

/* ---- Sampling kernel ----------------------------------------------------- */

static int _sample_one_primary(struct beam_workspace const *wb,
                               struct osh_rng *rng,
                               struct particle **part_out,
                               struct ray_v *ray_out) {
    struct beam_spot const *spot;
    double norm;
    double x, xp, y, yp;
    int rc;

    if (!wb || !rng || !part_out || !ray_out) {
        return OSH_EINVAL;
    }

    /* 1. Select spot */
    spot = _select_spot(wb);

    /* 2. Particle species */
    *part_out = spot->part;

    /* 3. Sample energy — stored in ray_out->p[3] */
    ray_out->p[3] = _sample_energy(rng, spot);

    /* 4. Sample transverse phase space (x, x') and (y, y') in PZALIGN */
    rc = _sample_position(rng, spot, &x, &y);
    if (rc != OSH_OK) {
        return rc;
    }
    _sample_angle_1d(rng, x, _shape_sigma_1d(spot, 0), spot->div[0], spot->cor[0], &xp);
    _sample_angle_1d(rng, y, _shape_sigma_1d(spot, 1), spot->div[1], spot->cor[1], &yp);

    ray_out->p[0] = x;
    ray_out->p[1] = y;
    ray_out->p[2] = 0.0; /* beam starts at z=0 in PZALIGN */

    /* Direction in PZALIGN: small-angle (x', y') plus main z component. */
    ray_out->v[0] = xp;
    ray_out->v[1] = yp;
    ray_out->v[2] = 1.0;

    ray_out->system = OSH_COORD_PZALIGN;

    /* 5. SAD correction */
    if (wb->shared.use_sad) {
        rc = _apply_sad(ray_out, spot, &wb->shared);
        if (rc != OSH_OK) {
            return rc;
        }
    }

    norm = sqrt(ray_out->v[0] * ray_out->v[0] + ray_out->v[1] * ray_out->v[1] + ray_out->v[2] * ray_out->v[2]);
    if (norm <= 0.0) {
        return OSH_ESTATE;
    }
    ray_out->v[0] /= norm;
    ray_out->v[1] /= norm;
    ray_out->v[2] /= norm;

    /* 6. Apply affine transform PZALIGN -> UNIVERSE */
    _apply_transform(ray_out, spot);

    return OSH_OK;
}

/* ---- Public entry points ------------------------------------------------- */

int osh_beam_new_primaries(
    struct beam_workspace const *wb, struct osh_rng *rng, size_t n, struct particle **part_out, struct ray_v *ray_out) {
    size_t i;
    int rc;

    if (!wb || !rng || !part_out || !ray_out) {
        return OSH_EINVAL;
    }

    /* Current implementation is intentionally simple: batch API first,
     * SoA/vectorized internals next. Keeping one looped kernel lets us evolve
     * toward SIMD/GPU without changing the caller-facing contract again. */
    for (i = 0; i < n; i++) {
        rc = _sample_one_primary(wb, rng, &part_out[i], &ray_out[i]);
        if (rc != OSH_OK) {
            return rc;
        }
    }

    return OSH_OK;
}

int osh_beam_new_primary(struct beam_workspace const *wb,
                         struct osh_rng *rng,
                         struct particle **part_out,
                         struct ray_v *ray_out) {
    return osh_beam_new_primaries(wb, rng, 1, part_out, ray_out);
}
