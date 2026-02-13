#ifndef _OSH_BEAM
#define _OSH_BEAM
#include <stdio.h>
#include <stdlib.h>

#include "common/osh_readline.h"
#include "particle/osh_particle.h"

#define OSH_BEAM_STRAGG_OFF 0
#define OSH_BEAM_STRAGG_GAUSS 1
#define OSH_BEAM_STRAGG_VAVILOV 2

#define OSH_BEAM_MSCAT_OFF 0
#define OSH_BEAM_MSCAT_GAUSS 1
#define OSH_BEAM_MSCAT_MOLIERE 2

#define OSH_BEAM_TMIN0 0.1 /* MeV or MeV/nucleon */

/* forward declarations */
struct ripple_filter;
struct parlev;
struct shfile;
struct ray_c;

/* a single beam spot */
struct beam_spot {
    /* position data is for nominal transport along Z */
    struct particle *part; /* pointer to what particle we have */
    double _tm[16];        /* matrix for beam rotation and translation */
    double p[3];           /* position [cm] relative to isocenter */
    double size[2];        /* beam size parameter (e.g. 1 sigma if gaussian, inner and outer radius if ciruclar) [cm] */
    double div[2];         /* beam divergence [rad] */
    double cov[2];         /* beam covariance [cm^2] */

    double t0;     /* initial energy [MeV] or [MeV/nucleon] */
    double tsigma; /* energy spread [MeV] or [MeV/nucleon] */
    double p0;     /* initial momentum [MeV/c] */
    double psigma; /* momentum spread [MeV/c] */

    double wt;             /* weight (1.0 = normal weight) (currently not fully by SHIELD-HIT core), leave it at 1.0 */
    unsigned int spot_id;  /* spot number */
    unsigned int layer_id; /* energy layer number */

    char shape;       /* shape of the beam, as listed in OSH_BEAM_SHAPE_* in osh_beamdef.h */
    char tsigma_type; /* energy spread distribution OSH_RANDOM_DIST_* (square, gaussian) */
};

struct beam_shared {
    double tcut[2]; /* upper and lower primary cutoff energy. */
    double pcut[2]; /* upper and lower primary cutoff momentum. */
    double sad[2];  /* scanning-magnet to isocenter distance [cm] */
    double focus;   /* focus point relative to source (k-distance) [cm] */
    /* theta and phi are by ISO 80000-2:2019 convention.
        see https://en.wikipedia.org/wiki/Spherical_coordinate_system

        direction axis : (theta, phi)
        x : (pi/2, 0)
        y : (pi/2,pi/2)
        z : (0,*)
     */
    double theta; /* beam theta relative to Z axis [rad]*/
    double phi;   /* beam phi realtive to Z axis [rad] */

    char use_pmax;   /* true if user specified momentum instead of energy  */
    char use_psigma; /* true if user specified momentum instead of energy */
    char use_div;    /* flag to skip divergence calculation, if not requested by user */
    char use_sad;    /* flag to skip SAD calculation, if not requested by user */
};

// /* --- Spot list --- */
// struct beam_spots {
//     struct beam_spot *spots; /* owned array */
//     struct beam_shared *shared; /* shared across all spots */
//     char *fname;  /* input file */
//     size_t nspots;
//     // double norm;
// };

/* a phase-space file */
struct beam_phsp {
    struct particle **part; /* list of pointers to particles */
    size_t len;             /* array length */
    double *p[3];           /* list of positions */
    double *d[3];           /* list of directions */
    double *wt;             /* list of weights */
    char *fname;            /* input file */
};

/* defining beam parameters */
struct beam_workspace {
    struct beam_phsp *phsp;     /* phase space files */
    struct beam_spot *spots;    /* placeholder for beam scanning data (struct is init to 0 if not given) */
    struct beam_shared *shared; /* shared parameters across all spots */
    struct ripple_filter *rifi; /* beam modulator placeholder (struct is init to 0 if not given) */
    struct parlev *parlev;      /* parameter settings from PARLEV() array found in /LEVERS/ common block */
    char *wdir;                 /* working directory */
    char *fname;                /* filename of beam input file */
    char *fname_spotlist;       /* filename of optional external spotlist */
    size_t nspots;              /* Number of spots */

    size_t nstat;   /* Number of requested primaries */
    size_t nsave;   /* saving step of primaries, 0 for disable */
    int rndseed;    /* random number seed as given in beam.dat */
    int rndoffset;  /* offset to rndseed as given by -N option */
    float deltae;   /* max. fraction of energy lost in a single step */
    float oln;      /* lower neutron energy cut */
    float demin;    /* lower moliere scatter cut */
    char straggl;   /* switch for energy straggling */
    char scatter;   /* switch for multiple scattering */
    char nuclear;   /* switch for nuclear reactions */
    char emtrans;   /* switch for EM-transport (unused) */
    char apcorr;    /* switch for alternative antiproton annihilation model */
    char beam_mode; /* switch for spots or PHSP*/
    char makeln;    /* switch for generation of neutron output data file */
    char neutrfast; /* switch for fast neutron transport */
};

/**
 * @brief Allocate, parse, and fully initialize a beam workspace.
 *
 * This is the one-shot, library-style entry point. It:
 *  - allocates and zero-initializes a new ::beam_workspace,
 *  - applies sane defaults,
 *  - parses @p filename (using @p wdir as working directory if non-NULL),
 *  - builds the spot list or phase-space (depending on the file contents),
 *  - wires shared, run-wide parameters (spotlist->shared = &wb->shared),
 *  - computes per-spot transform matrices (spot->_tm[16]),
 *  - validates invariants (e.g., non-empty spotlist/PHSP, sane cutoffs),
 *  - and returns a ready-to-use workspace via @p wb_out.
 *
 * On success, @p *wb_out owns all allocated resources and must be released with
 * ::osh_beam_workspace_free(). On failure, no allocation is leaked and @p *wb_out
 * is left unmodified.
 *
 * ### Inputs and ownership
 * - @p filename: path to a configuration file. The content defines either a spot list
 *   (beam_mode=spots) or a phase-space source (beam_mode=phsp). The string is not retained.
 * - @p wdir: optional working directory for relative includes/paths; copied internally
 *   if provided; may be NULL.
 * - @p wb_out: non-NULL pointer to receive the allocated workspace pointer.
 *
 * ### Outputs and guarantees
 * - @p *wb_out is non-NULL on success and points to a fully-initialized ::beam_workspace.
 * - If spots are configured:
 *     - wb->spotlist is non-NULL, wb->spotlist->nspots > 0,
 *     - wb->spotlist->shared == &wb->shared (wired),
 *     - every spot has its per-spot transform filled in spot->_tm[16].
 * - If PHSP is configured:
 *     - wb->phsp is non-NULL, wb->phsp->len > 0.
 * - Units: energies are MeV (or MeV/u if per-nucleon), momenta are MeV/c (or MeV/c/u).
 *   The exact interpretation of t0/p0 is preserved from the input (see struct comments).
 *
 * @return error code as in osh_err.h
 *
 * @note Thread-safety: this function is not thread-safe. Use per-thread workspaces.
 *
 * @param filename  Path to the configuration file (not retained).
 * @param wdir      Optional working directory; may be NULL (copied if non-NULL).
 * @param wb_out    Output: on success, set to a newly allocated ::beam_workspace.
 *
 * @return 0 on success; negative OSH_E* code on failure (see osh_err.h).
 *
 * @author Niels Bassler
 */
int osh_beam_setup(char const *filename, char const *wdir, struct beam_workspace **wb);

/**
 * @brief Release a beam workspace and all owned resources.
 *
 * Frees all memory owned by @p wb, including nested allocations (spot arrays,
 * PHSP arrays, strings, etc.). Safe to call with NULL.
 *
 * @param wb  Workspace returned by ::osh_beam_setup(), or NULL.
 * @return 0 (OSH_OK). Always succeeds.
 *
 * @author Niels Bassler
 */
int osh_beam_workspace_free(struct beam_workspace *wb);
// int osh_beam_get_primary(struct beam_workspace *wb, struct particle **part, struct ray_c *ray);  TODO later

/**
 * @brief Print a concise summary of the workspace to stdout.
 *
 * Intended for debugging/logging. Does not modify @p wb.
 *
 * @param wb  Workspace to print; ignored if NULL.
 *
 * @author Niels Bassler
 */
void osh_beam_print(struct beam_workspace const *wb);

/**
 * @brief Print a concise summary of a single beam spot to stdout.
 *
 * Shows position, size/divergence, t0/p0 with spreads, IDs, and shape.
 * Does not modify @p spot.
 *
 * @param spot  Spot to print; ignored if NULL.
 *
 * @author Niels Bassler
 */
void osh_beam_print_spot(struct beam_spot const *spot);

#endif /* !_OSH_BEAM */