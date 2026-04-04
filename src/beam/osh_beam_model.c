#include "beam/osh_beam_model.h"

#include <math.h>

#include "common/osh_const.h"
#include "common/osh_coord.h"
#include "common/osh_interpolate.h"
#include "common/osh_rc.h"
#include "common/osh_vect.h"

/* ---- Step 1: spot selection ---------------------------------------------- */

/* Return a pointer to the spot to use for this history.
 * Single-spot and phase-space modes always return spots[0].
 * SOBP mode draws a spot from the list weighted by spot->wt. */
static struct beam_spot const *_select_spot(struct beam_workspace const *wb, struct osh_rng *rng) {
    double w;
    long int idx;

    if (!wb || !wb->spots || wb->nspots == 0) {
        return NULL;
    }
    if (wb->nspots == 1) {
        return &wb->spots[0];
    }
    if (!rng || !wb->cum_wt || wb->wt_sum <= 0.0) {
        return NULL;
    }

    w = osh_rng_double(rng) * wb->wt_sum;
    idx = osh_binary_search_upper_d(w, wb->cum_wt, wb->nspots);
    if (idx < 0 || idx >= (long int) wb->nspots) {
        return NULL;
    }
    return &wb->spots[idx];
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

static int _sample_position_xy(struct osh_rng *rng, struct beam_spot const *spot, double pos[2]) {
    double rin;
    double rout;
    double r2;
    double phi;

    switch (spot->shape) {
    case OSH_BEAM_SHAPE_PENCIL:
        pos[0] = 0.0;
        pos[1] = 0.0;
        return OSH_OK;

    case OSH_BEAM_SHAPE_GAUSSIAN:
        pos[0] = (spot->size[0] > 0.0) ? spot->size[0] * osh_rng_gauss01(rng) : 0.0;
        pos[1] = (spot->size[1] > 0.0) ? spot->size[1] * osh_rng_gauss01(rng) : 0.0;
        return OSH_OK;

    case OSH_BEAM_SHAPE_SQUARE:
        pos[0] = _sample_uniform_symmetric(rng, spot->size[0]);
        pos[1] = _sample_uniform_symmetric(rng, spot->size[1]);
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
        pos[0] = rout * cos(phi);
        pos[1] = rout * sin(phi);
        return OSH_OK;

    default:
        return OSH_EINVAL;
    }
}

/* Sample the complete transverse phase space in beam-local coordinates.
 *
 * Output:
 *   pos[0], pos[1] = sampled transverse offsets x,y at the beam start plane
 *   ang[0], ang[1] = sampled small-angle slopes x'=dx/dz, y'=dy/dz
 *
 * The current beam model treats the two transverse planes independently:
 *   - size[0], div[0], cor[0] describe the x/x' plane
 *   - size[1], div[1], cor[1] describe the y/y' plane
 *   - there are no cross-plane x-y coupling terms yet
 *
 * Gaussian spots sample position and angle together in each plane so the
 * requested position-angle correlation is exact. Non-Gaussian shapes first
 * sample the transverse position from the requested profile, then derive the
 * angle using the RMS width of that profile in the same correlation formula. */
static int _sample_phase_space_xy(struct osh_rng *rng, struct beam_spot const *spot, double pos[2], double ang[2]) {
    double sigma_pos[2];
    double rho[2];
    double rho_orth[2];
    double z_ang[2];
    int axis;
    int rc;

    if (spot->shape == OSH_BEAM_SHAPE_GAUSSIAN) {
        for (axis = 0; axis < 2; axis++) {
            /* In the Gaussian case, position and angle belong to one joint
             * 2x2 phase-space model in each plane, so sample them together. */
            sigma_pos[axis] = _shape_sigma_1d(spot, axis);
            pos[axis] = (sigma_pos[axis] > 0.0) ? sigma_pos[axis] * osh_rng_gauss01(rng) : 0.0;

            if (spot->div[axis] <= 0.0) {
                ang[axis] = 0.0;
                continue;
            }

            if (sigma_pos[axis] <= 0.0) {
                ang[axis] = spot->div[axis] * osh_rng_gauss01(rng);
                continue;
            }

            rho[axis] = spot->cor[axis];
            if (rho[axis] > 1.0) {
                rho[axis] = 1.0;
            } else if (rho[axis] < -1.0) {
                rho[axis] = -1.0;
            }

            rho_orth[axis] = sqrt(fmax(0.0, 1.0 - rho[axis] * rho[axis]));
            z_ang[axis] = osh_rng_gauss01(rng);
            /* Closed-form correlated Gaussian:
             * x' = sigma_x' * (rho * z_pos + sqrt(1-rho^2) * z_ang)
             * with z_pos = x / sigma_x. */
            ang[axis] = spot->div[axis] * (rho[axis] * (pos[axis] / sigma_pos[axis]) + rho_orth[axis] * z_ang[axis]);
        }
        return OSH_OK;
    }

    /* For square/circular/pencil beams the spatial profile is not Gaussian.
     * Sample the position from the requested shape first, then map that
     * sampled position to an angle using the RMS width of the profile. */
    rc = _sample_position_xy(rng, spot, pos);
    if (rc != OSH_OK) {
        return rc;
    }

    for (axis = 0; axis < 2; axis++) {
        sigma_pos[axis] = _shape_sigma_1d(spot, axis);
        if (spot->div[axis] <= 0.0) {
            ang[axis] = 0.0;
            continue;
        }

        if (sigma_pos[axis] <= 0.0) {
            ang[axis] = spot->div[axis] * osh_rng_gauss01(rng);
            continue;
        }

        rho[axis] = spot->cor[axis];
        if (rho[axis] > 1.0) {
            rho[axis] = 1.0;
        } else if (rho[axis] < -1.0) {
            rho[axis] = -1.0;
        }

        rho_orth[axis] = sqrt(fmax(0.0, 1.0 - rho[axis] * rho[axis]));
        z_ang[axis] = osh_rng_gauss01(rng);
        /* Same correlation formula as above, but with the actual sampled
         * position and the RMS width implied by the chosen beam shape. */
        ang[axis] = spot->div[axis] * (rho[axis] * (pos[axis] / sigma_pos[axis]) + rho_orth[axis] * z_ang[axis]);
    }

    return OSH_OK;
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
    double pos[2];
    double ang[2];
    double norm;
    int rc;

    if (!wb || !rng || !part_out || !ray_out) {
        return OSH_EINVAL;
    }

    /* 1. Select spot */
    spot = _select_spot(wb, rng);
    if (!spot) {
        return OSH_ESTATE;
    }

    /* 2. Particle species */
    *part_out = spot->part;

    /* 3. Sample energy — stored in ray_out->p[3] */
    ray_out->p[3] = _sample_energy(rng, spot);

    /* 4. Sample transverse phase space (x, y, x', y') in PZALIGN */
    rc = _sample_phase_space_xy(rng, spot, pos, ang);
    if (rc != OSH_OK) {
        return rc;
    }

    ray_out->p[0] = pos[0];
    ray_out->p[1] = pos[1];
    ray_out->p[2] = 0.0; /* beam starts at z=0 in PZALIGN */

    /* Direction in PZALIGN: small-angle (x', y') plus main z component. */
    ray_out->v[0] = ang[0];
    ray_out->v[1] = ang[1];
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
