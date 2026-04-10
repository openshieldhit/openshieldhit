#include "beam/osh_beam_model.h"

#include <math.h>

#include "common/osh_const.h"
#include "common/osh_coord.h"
#include "common/osh_interpolate.h"
#include "common/osh_rc.h"
#include "common/osh_vect.h"

/* ---- Step 1: spot selection ---------------------------------------------- */

/**
 * @brief Return a pointer to the beam spot to use for this primary history.
 *
 * @details
 * Single-spot beams always return spots[0]. SOBP mode draws from the spot
 * list using inverse-CDF sampling on the cumulative weight array cum_wt:
 * a uniform deviate scaled to [0, wt_sum] is located with an upper-bound
 * binary search, giving each spot a selection probability proportional to
 * its weight.
 *
 * @param[in] wb   Beam workspace with a fully initialised spot list and
 *                 cumulative weight array.
 * @param[in] rng  Random-number generator; only consumed in SOBP mode.
 *
 * @returns Pointer into wb->spots, or NULL on invalid input.
 */
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

/**
 * @brief Sample total kinetic energy [MeV] for one primary.
 *
 * @details
 * Draws from a Gaussian N(t0, tsigma²) and clamps the result to zero to
 * avoid unphysical negative energies. When tsigma <= 0 the beam is
 * mono-energetic and t0 is returned directly.
 *
 * @param[in] rng   Random-number generator.
 * @param[in] spot  Beam spot carrying t0 [MeV] and tsigma [MeV].
 *
 * @returns Sampled kinetic energy in MeV, >= 0.
 */
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

/**
 * @brief Sample from a symmetric uniform distribution U(-half_width, +half_width).
 *
 * @param[in] rng        Random-number generator.
 * @param[in] half_width Half-width of the interval.
 *
 * @returns A value in (-half_width, +half_width).
 */
static double _sample_uniform_symmetric(struct osh_rng *rng, double half_width) {
    return (2.0 * osh_rng_double(rng) - 1.0) * half_width;
}

/**
 * @brief Return the RMS width of the beam profile along one transverse axis.
 *
 * @details
 * Used in the position-angle correlation formula as sigma_pos, i.e. the
 * denominator in  ang = div * (rho * pos/sigma_pos + sqrt(1-rho²) * z).
 * For Gaussian beams sigma_pos is also used to draw the position itself.
 *
 * Derivations per shape:
 *
 *   PENCIL   — zero width by definition; sigma = 0.
 *
 *   GAUSSIAN — size[axis] is already the Gaussian sigma; sigma = size[axis].
 *
 *   SQUARE   — uniform on [-w, +w] where w = size[axis].
 *              Var[x] = (1/2w) * integral_{-w}^{w} x² dx = w²/3,
 *              so sigma = w / sqrt(3).
 *
 *   CIRCULAR — uniform annulus with inner radius r_in, outer radius r_out.
 *              E[x²] = (1 / pi*(r_out²-r_in²))
 *                      * integral_0^{2pi} cos²(theta) d(theta)
 *                      * integral_{r_in}^{r_out} r³ dr
 *                    = (r_out⁴ - r_in⁴) / (4*(r_out²-r_in²))
 *                    = (r_out² + r_in²) / 4,
 *              so sigma = 0.5 * sqrt(r_in² + r_out²).
 *              The parser always sets size[1]=0 (solid disk), giving
 *              sigma = r_out / 2.
 *
 * @param[in] spot  Beam spot carrying shape and size[].
 * @param[in] axis  Transverse axis index: 0 = x, 1 = y.
 *
 * @returns RMS sigma >= 0, or -1.0 for an unrecognised shape.
 */
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

/**
 * @brief Sample the transverse position (x, y) from the beam spot's spatial profile.
 *
 * @details
 * Only called for non-Gaussian shapes; the Gaussian case is handled inline
 * in _sample_phase_space_xy to allow joint position-angle sampling.
 *
 * PENCIL   — pos = (0, 0).
 * SQUARE   — two independent uniform draws on [-size[0], +size[0]] and
 *            [-size[1], +size[1]].
 * CIRCULAR — polar sampling: phi uniform on [0, 2*pi), r² uniform on
 *            [r_in², r_out²], giving a uniform density over the annular area.
 *
 * @param[in]  rng   Random-number generator.
 * @param[in]  spot  Beam spot with shape and size[].
 * @param[out] pos   Sampled position {x, y} in cm (beam-local PZALIGN frame).
 *
 * @returns OSH_OK on success, OSH_EINVAL for an invalid or unrecognised shape.
 */
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

/**
 * @brief Sample the complete transverse phase space in beam-local coordinates.
 *
 * @details
 * The two transverse planes (x/x' and y/y') are treated independently;
 * there are no cross-plane coupling terms.
 *
 * For each plane the angle is drawn from the correlated Gaussian model:
 *
 *   ang = sigma_ang * (rho * z_pos + sqrt(1 - rho²) * z_ang)
 *
 * where z_pos = pos / sigma_pos and z_ang are independent N(0,1) draws, rho
 * is the position-angle correlation coefficient (spot->cor[axis]), and
 * sigma_ang = spot->div[axis].  This gives Cor(pos, ang) = rho exactly,
 * regardless of the spatial profile shape, because E[pos * ang] / (sigma_pos
 * * sigma_ang) = rho * E[z_pos²] = rho.
 *
 * Gaussian spots: position and angle are sampled together in one step so
 * both marginals are exactly Gaussian and the joint distribution is a proper
 * bivariate Gaussian (emittance ellipse).
 *
 * Non-Gaussian spots (SQUARE, CIRCULAR): position is sampled first from the
 * requested spatial profile, then the angle is derived via the formula above
 * using _shape_sigma_1d() for sigma_pos.  The angular marginal is still
 * Gaussian; the joint distribution is not bivariate Gaussian but preserves
 * the correct RMS divergence and correlation coefficient.
 *
 * @param[in]  rng   Random-number generator.
 * @param[in]  spot  Beam spot with shape, size[], div[], and cor[].
 * @param[out] pos   Sampled transverse offsets {x, y} [cm] in PZALIGN frame.
 * @param[out] ang   Sampled small-angle slopes {x'=dx/dz, y'=dy/dz} [rad].
 *
 * @returns OSH_OK on success, OSH_EINVAL for an unrecognised shape.
 */
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

/**
 * @brief Apply source-axis-distance (SAD) fan-out correction in PZALIGN.
 *
 * @details
 * SAD models a scanning magnet whose virtual source lies a distance sad[axis]
 * upstream of isocenter.  The fan-out angle for each transverse plane is:
 *
 *   ray->v[axis] += (spot->p[axis] + ray->p[axis]) / (spot->p[2] + ray->p[2] + sad[axis])
 *
 * where the numerator is the full transverse offset at the start plane and
 * the denominator is the downstream distance from the virtual source.
 * BEAMPOS is already folded into spot->p by the parser, so this function
 * updates only the direction vector, not the start position.
 *
 * @param[in,out] ray  Ray in PZALIGN; direction updated in place.
 * @param[in]     spot Current beam spot (provides the nominal start position).
 * @param[in]     sh   Shared beam parameters carrying sad[0] and sad[1] [cm].
 *
 * @returns OSH_OK on success, OSH_ESTATE if the virtual-source distance is
 *          non-positive (degenerate geometry).
 */
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

/**
 * @brief Apply the precomputed affine transform PZALIGN -> UNIVERSE.
 *
 * @details
 * Executes the standard rigid-body transform stored in spot->_tm[16]:
 *
 *   p_universe = R * p_local + t
 *   v_universe = R * v_local
 *
 * This intentionally does NOT use osh_ray_v_transform(), which follows the
 * legacy SHIELD-HIT/GEMCA sign convention for the translation component.
 * The matrix spot->_tm is built by the post-parse step and already folds in
 * the beam direction (theta, phi) and the BEAMPOS offset.
 *
 * @param[in,out] ray   Ray in PZALIGN on entry; converted to UNIVERSE in place.
 * @param[in]     spot  Provides the 4x4 row-major affine matrix _tm[16],
 *                      as consumed by osh_vect_trans_point_affine() and
 *                      osh_vect_trans_vector_affine().
 */
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

/**
 * @brief Sample one complete primary particle in UNIVERSE coordinates.
 *
 * @details
 * Executes the full five-step sampling pipeline:
 *   1. Select beam spot (weighted draw for SOBP, spots[0] otherwise).
 *   2. Sample kinetic energy around spot->t0 with spread spot->tsigma.
 *   3. Sample transverse phase space (x, y, x', y') in PZALIGN.
 *   4. Apply SAD fan-out correction if shared.use_sad is set.
 *   5. Normalise the direction vector and apply the affine transform to UNIVERSE.
 *
 * @param[in]  wb        Fully initialised beam workspace.
 * @param[in]  rng       Random-number generator.
 * @param[out] part_out  Receives a pointer to the sampled particle species
 *                       (owned by wb; caller must not free).
 * @param[out] ray_out   Receives the sampled ray in OSH_COORD_UNIVERSE;
 *                       ray_out->p[3] holds total kinetic energy [MeV].
 *
 * @returns OSH_OK on success, OSH_E* on error.
 */
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
