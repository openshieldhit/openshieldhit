#include "beam/osh_beam_model.h"

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
 * TODO: implement Gaussian draw using spot->tsigma and spot->tsigma_type. */
static double _sample_energy(struct beam_spot const *spot) {
    /* TODO: draw from distribution around spot->t0 with spread spot->tsigma */
    return spot->t0;
}

/* ---- Step 3: transverse phase-space sampling ----------------------------- */

/* Sample (position, angle) in one transverse plane using the Fermi-Eyges
 * Cholesky form.  sigma_x and sigma_xp are the RMS widths; rho is the
 * correlation coefficient.  u1 and u2 are independent N(0,1) samples.
 *
 *   x  = sigma_x  * u1
 *   x' = sigma_xp * (rho * u1 + sqrt(1 - rho^2) * u2)
 *
 * TODO: plug in actual random number generator; handle SQUARE/CIRCULAR shapes. */
static void
_sample_phasespace_1d(struct osh_rng *rng, double sigma_x, double sigma_xp, double rho, double *x_out, double *xp_out) {
    /* TODO: replace with actual RNG calls */
    double u1 = 0.0;
    double u2 = 0.0;
    (void) rng;
    (void) rho;

    *x_out = sigma_x * u1;
    *xp_out = sigma_xp * u2; /* TODO: full Cholesky: rho*u1 + sqrt(1-rho^2)*u2 */
}

/* ---- Step 4: SAD correction ---------------------------------------------- */

/* Apply source-axis-distance correction in PZALIGN.
 * TODO: implement fan-out shift and angle update. */
static void _apply_sad(struct ray_v *ray, struct beam_spot const *spot, struct beam_shared const *sh) {
    (void) ray;
    (void) spot;
    (void) sh;
    /* TODO */
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
    double x, xp, y, yp;

    if (!wb || !part_out || !ray_out) {
        return OSH_EINVAL;
    }

    /* 1. Select spot */
    spot = _select_spot(wb);

    /* 2. Particle species */
    *part_out = spot->part;

    /* 3. Sample energy — stored in ray_out->p[3] */
    ray_out->p[3] = _sample_energy(spot);

    /* 4. Sample transverse phase space (x, x') and (y, y') in PZALIGN */
    _sample_phasespace_1d(rng, spot->size[0], spot->div[0], spot->cov[0], &x, &xp);
    _sample_phasespace_1d(rng, spot->size[1], spot->div[1], spot->cov[1], &y, &yp);

    ray_out->p[0] = x;
    ray_out->p[1] = y;
    ray_out->p[2] = 0.0; /* beam starts at z=0 in PZALIGN */

    /* Direction in PZALIGN: small-angle (x', y') plus main z component */
    ray_out->v[0] = xp;
    ray_out->v[1] = yp;
    ray_out->v[2] = 1.0; /* TODO: normalise after setting xp, yp */

    ray_out->system = OSH_COORD_PZALIGN;

    /* 5. SAD correction */
    if (wb->shared.use_sad) {
        _apply_sad(ray_out, spot, &wb->shared);
    }

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
