#include "beam/osh_beam_model.h"

#include <math.h>

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
 * Current implementation follows the beam-model path used in export_mcpl.py:
 * Gaussian sampling around t0 with sigma tsigma, clamped at zero. */
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

/* Sample Gaussian (position, angle) in one transverse plane using the
 * Fermi-Eyges Cholesky form. sigma_x and sigma_xp are the RMS widths; rho is
 * the correlation coefficient. u1 and u2 are independent N(0,1) samples.
 *
 *   x  = sigma_x  * u1
 *   x' = sigma_xp * (rho * u1 + sqrt(1 - rho^2) * u2)
 *
 * The Python reference stores rho, not covariance, and computes the equivalent
 * 2x2 Cholesky factors once per beam layer. For the current single-spot path
 * we use the equivalent closed form directly. */
static void _sample_phasespace_gaussian_1d(
    struct osh_rng *rng, double sigma_x, double sigma_xp, double rho, double *x_out, double *xp_out) {
    double u1;
    double u2;
    double rho_clamped;
    double rho_orth;

    if (sigma_x <= 0.0) {
        *x_out = 0.0;
    } else {
        u1 = osh_rng_gauss01(rng);
        *x_out = sigma_x * u1;
        if (sigma_xp <= 0.0) {
            *xp_out = 0.0;
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
        *xp_out = sigma_xp * (rho_clamped * u1 + rho_orth * u2);
        return;
    }

    if (sigma_xp <= 0.0) {
        *xp_out = 0.0;
    } else {
        *xp_out = sigma_xp * osh_rng_gauss01(rng);
    }
}

static int
_sample_phasespace_1d(struct osh_rng *rng, struct beam_spot const *spot, int axis, double *x_out, double *xp_out) {
    switch (spot->shape) {
    case OSH_BEAM_SHAPE_PENCIL:
        *x_out = 0.0;
        *xp_out = 0.0;
        return OSH_OK;

    case OSH_BEAM_SHAPE_GAUSSIAN:
        _sample_phasespace_gaussian_1d(rng, spot->size[axis], spot->div[axis], spot->cor[axis], x_out, xp_out);
        return OSH_OK;

    case OSH_BEAM_SHAPE_SQUARE:
    case OSH_BEAM_SHAPE_CIRCULAR:
        return OSH_ENOTSUP;

    default:
        return OSH_EINVAL;
    }
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
    rc = _sample_phasespace_1d(rng, spot, 0, &x, &xp);
    if (rc != OSH_OK) {
        return rc;
    }
    rc = _sample_phasespace_1d(rng, spot, 1, &y, &yp);
    if (rc != OSH_OK) {
        return rc;
    }

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
