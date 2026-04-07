#ifndef OSH_BEAM_H
#define OSH_BEAM_H

#include <stddef.h>
#include <stdint.h>

#include "beam/osh_beamdef.h"
#include "common/osh_logger.h"
#include "particle/osh_particle.h"

/* ---- Forward declarations ------------------------------------------------ */

struct ripple_filter; /* beam range modulator — defined in osh_beam_rifi.h */
struct parlev;        /* LEVERS parameter block — defined in osh_beam_parlev.h */

/* ---- beam_spot ------------------------------------------------------------ */

/* A single beam spot.
 *
 * After osh_beam_setup_from_path(), all energies are TOTAL kinetic energy in
 * MeV — NOT MeV/nucleon, NOT specific energy. User-facing beam inputs such as
 * TMAX0 use MeV/nucleon for ions and are converted during post-parse once the
 * primary particle is known. The transport engine divides by part->a when it
 * needs specific energy for stopping-power table lookup. Storing total values
 * here avoids per-primary multiplications in the inner transport loop.
 *
 * All momenta are total momentum in MeV/c — NOT per nucleon.
 *
 * Fields prefixed with '_' are derived quantities computed once by
 * osh_beam_setup_from_path() after parsing. The parser and transport loop must
 * not write to them. */
struct beam_spot {
    /* Derived at init — do not write after osh_beam_setup_from_path() returns */
    double _tm[16]; /* 4x4 standard affine matrix mapping beam-local PZALIGN
                     * coordinates to UNIVERSE for sampling:
                     *   p_u = R * p_l + R * spot->p
                     * where spot->p is BEAMPOS in the local beam frame before
                     * gantry/table rotation. The 3x3 block rotates local beam
                     * basis vectors into UNIVERSE; the last column stores the
                     * translated world-space beam entry point. */

    /* Lateral profile */
    double p[3];    /* BEAMPOS [cm] in beam-local coordinates before gantry /
                     * table rotation; relative to isocenter */
    double size[2]; /* lateral size [cm]: 1-sigma (Gaussian), half-width (square),
                       inner/outer radius (circular) */
    double div[2];  /* angular spread [rad] */
    double cor[2];  /* position-angle correlation coefficient rho in [-1,1] */

    /* Energy / momentum — total quantities, NOT per nucleon */
    double t0;     /* total kinetic energy [MeV] */
    double tsigma; /* total kinetic energy spread [MeV] */
    double p0;     /* total momentum [MeV/c]; derived from t0 at init if not given explicitly */
    double psigma; /* total momentum spread [MeV/c]; derived from tsigma at init if not given */

    /* Sampling metadata */
    double wt; /* statistical weight (1.0 = normal) */
    unsigned int spot_id;
    unsigned int layer_id;

    struct particle *part;   /* particle species — pointer into a shared table, NOT owned here */
    char shape;              /* lateral profile shape: OSH_BEAM_SHAPE_* */
    char tsigma_type;        /* energy spread distribution: OSH_RANDOM_DIST_* */
    char t0_per_nucleon;     /* parser input marker; post-parse converts t0 to total MeV */
    char tsigma_per_nucleon; /* parser input marker; post-parse converts tsigma to total MeV */
};

/* ---- beam_shared ---------------------------------------------------------- */

/* Beam source geometry shared across all spots in a beam file.
 *
 * This struct covers beam direction, source-to-axis geometry, divergence flags,
 * and derived per-run maxima used for stopping-power table sizing.
 *
 * Transport control parameters (cutoffs, step limits, physics switches) are
 * NOT here — they live in beam_workspace.
 *
 * emax and pmax are computed once by osh_beam_setup_from_path() as the maximum
 * over all spots. Units: total MeV and total MeV/c respectively. */
struct beam_shared {
    /* Derived at init — do not write after osh_beam_setup_from_path() returns */
    double emax; /* max t0 over all spots [MeV] — for stopping-power table cap */
    double pmax; /* max p0 over all spots [MeV/c] */

    double sad[2]; /* source-to-isocenter / focal distance [cm]: [0]=x, [1]=y.
                    * This is a positive machine-geometry distance, not a
                    * signed coordinate. The corresponding virtual source is
                    * always upstream of isocenter along the local beam axis.
                    * SAD is therefore defined relative to isocenter, NOT
                    * relative to BEAMPOS or the current beam start plane. */
    double focus;  /* focus point along beam axis relative to source [cm] */

    /* Beam direction by ISO 80000-2:2019 spherical convention:
     *   theta=0            -> along +Z  (any phi)
     *   theta=pi/2, phi=0  -> along +X
     *   theta=pi/2, phi=pi/2 -> along +Y */
    double theta; /* polar angle from +Z [rad], range [0, pi] */
    double phi;   /* azimuthal angle from +X [rad], range [0, 2*pi) */

    char use_div; /* 1 if divergence sampling is active */
    char use_sad; /* 1 if SAD correction is active */
};

/* ---- beam_phsp ------------------------------------------------------------ */

/* Phase-space particle source, following the MCPL file convention:
 *   - direction vectors are normalised unit vectors
 *   - energy is total kinetic energy in MeV (NOT per nucleon)
 *   - particle species is stored as PDG particle code per entry
 *
 * Positions and directions use structure-of-arrays layout (p[3], d[3]) so each
 * coordinate axis is a contiguous double array, enabling SIMD vectorisation.
 *
 * The particle struct pointer is not stored per entry to save memory. A
 * one-entry cache (_cached_pdg / _cached_part) avoids repeated table lookups:
 * homogeneous files (single species) pay the lookup cost exactly once; mixed
 * files pay it on each species transition. */
struct beam_phsp {
    double *p[3]; /* position [cm]: p[0]=x, p[1]=y, p[2]=z — each a contiguous array */
    double *d[3]; /* normalised direction: d[0]=ux, d[1]=uy, d[2]=uz */
    double *e;    /* total kinetic energy per entry [MeV] — NOT per nucleon */
    double *wt;   /* statistical weight per entry */
    int32_t *pdg; /* PDG particle code per entry */
    size_t len;   /* number of entries */
    char *fname;  /* source file path (owned) */

    /* Sampling cache — avoids repeated particle table lookup */
    int32_t _cached_pdg;
    struct particle *_cached_part;
};

/* ---- beam_workspace ------------------------------------------------------- */

/* Top-level beam workspace.
 *
 * Holds the beam source (spots or phase-space), beam geometry (via shared),
 * and simulation-wide transport defaults parsed from beam.dat.
 *
 * nstat may be overridden after osh_beam_setup_from_path() by the caller if an
 * explicit CLI value was given (check has_nstat / cfg.nstat in openshieldhit_run).
 *
 * Transport cutoffs are global defaults for the run. Per-medium overrides are
 * stored in struct zone (osh_gemca2.h) and take precedence at simulation time. */
struct beam_workspace {
    /* --- Beam source --- */
    struct beam_phsp *phsp;     /* phase-space source; NULL if spots mode */
    struct beam_spot *spots;    /* spot array; NULL if phsp mode */
    struct particle primary;    /* resolved primary species shared by all spots */
    struct beam_shared shared;  /* beam geometry — embedded by value, not a pointer.
                                 * Each workspace is an independent run; there is no sharing
                                 * scenario that would justify a separate allocation. Embedding
                                 * avoids an extra malloc/free pair and keeps shared fields in
                                 * the same cache line as the rest of the workspace scalars. */
    struct ripple_filter *rifi; /* range modulator; NULL if not used */
    struct parlev *parlev;      /* PARLEV lever settings; NULL if not used */
    char *wdir;                 /* working directory for relative paths (owned) */
    char *fname;                /* beam input file path (owned) */
    char *fname_spotlist;       /* external spot list file (owned); NULL if not used */
    size_t nspots;              /* number of entries in spots[] */
    double *cum_wt;             /* cumulative spot weights for weighted SOBP selection */
    double wt_sum;              /* total spot weight; equals cum_wt[nspots-1] when nspots > 0 */
    char has_primary;           /* 1 once PRIMARY resolved successfully */

    /* --- Run control (may be overridden by CLI after setup) --- */
    size_t nstat;  /* number of requested primary histories */
    size_t nsave;  /* history save interval (0 = disabled) */
    int rndseed;   /* random number seed */
    int rndoffset; /* seed offset applied on top of rndseed */

    /* --- Transport cutoffs [global defaults; per-medium overrides in geometry] --- */
    float tcut;   /* lower primary ion energy cutoff [MeV/nucleon] */
    float pcut;   /* lower primary momentum cutoff [MeV/c] */
    float ncut;   /* lower neutron energy cutoff [MeV] */
    float deltae; /* max fractional energy loss per step */
    float demin;  /* lower Moliere multiple-scattering cutoff [MeV/nucleon] */

    /* --- Physics switches --- */
    char straggl;   /* energy straggling model: OSH_BEAM_STRAGG_* */
    char scatter;   /* multiple scattering model: OSH_BEAM_MSCAT_* */
    char nuclear;   /* nuclear reactions: 0=off, 1=on */
    char emtrans;   /* EM transport (unused, reserved for future use) */
    char apcorr;    /* alternative antiproton annihilation model */
    char beam_mode; /* source type: OSH_BEAM_MODE_* */
    char makeln;    /* generate neutron output file */
    char neutrfast; /* enable fast neutron transport */
};

/* ---- API ------------------------------------------------------------------ */

/**
 * @brief Allocate, parse, and fully initialise a beam workspace from a file.
 *
 * @details
 * The working directory for resolving relative paths inside the file (spot
 * lists, USEBMOD, USEPARLEV, etc.) is derived from dirname(path) automatically.
 * The sequence is: open file → parse raw keys → load external spot list (if
 * any) → validate → post-parse (derive energies, build affine matrices,
 * accumulate cumulative weights).
 *
 * On success *wb_out owns all allocated resources and must be released with
 * osh_beam_workspace_free(). On failure *wb_out is left unmodified and no
 * memory is leaked.
 *
 * Future constructor variants follow the same signature pattern:
 *   osh_beam_setup_from_pipe(fd,   lg, wb_out)
 *   osh_beam_setup_from_json(json, lg, wb_out)
 *
 * @param[in]  path    Path to the beam input file (not retained after return).
 * @param[in]  lg      Logger for diagnostics; NULL uses the global default.
 * @param[out] wb_out  Receives the allocated workspace pointer on success.
 *
 * @returns OSH_OK on success, OSH_E* on failure.
 */
int osh_beam_setup_from_path(char const *path, struct osh_logger *lg, struct beam_workspace **wb_out);

/**
 * @brief Release a beam workspace and all owned resources.
 *
 * @details
 * Frees spots, cum_wt, phsp, rifi, parlev, wdir, fname, and fname_spotlist.
 * The embedded beam_shared is not heap-allocated and needs no separate free.
 * Safe to call with NULL.
 *
 * @param[in] wb  Workspace to release; may be NULL.
 *
 * @returns OSH_OK always.
 */
int osh_beam_workspace_free(struct beam_workspace *wb);

/**
 * @brief Release a phase-space source and all its owned arrays.
 *
 * @details
 * Frees the six SoA arrays (p[0..2], d[0..2]), e, wt, pdg, and fname,
 * then frees the struct itself. Safe to call with NULL.
 *
 * @param[in] phsp  Phase-space source to release; may be NULL.
 *
 * @returns OSH_OK always.
 */
int osh_beam_phsp_free(struct beam_phsp *phsp);

/* TODO: osh_beam_get_primary() — sample next primary ray from workspace */

/**
 * @brief Print a concise summary of the beam workspace to the info logger.
 *
 * @param[in] wb  Workspace to print; silently ignored when NULL.
 */
void osh_beam_print(struct beam_workspace const *wb);

/**
 * @brief Print a concise summary of one beam spot to the info logger.
 *
 * @param[in] spot  Spot to print; silently ignored when NULL.
 */
void osh_beam_print_spot(struct beam_spot const *spot);

#endif /* OSH_BEAM_H */
