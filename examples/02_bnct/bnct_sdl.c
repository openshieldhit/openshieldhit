/* simple temporary test file for playing with the new scorer */
// gcc SDLtest.c -o dsl -lSDL2 -Wall
#include <SDL.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "common/osh_coord.h"
#include "common/osh_logger.h"
#include "common/osh_vect.h"
#include "common/osh_version.h"
#include "gemca/osh_gemca2.h"
#include "random/osh_rng.h"
#include "transport/osh_transport.h"

#define WINDOW_WIDTH 800
#define WINDOW_HEIGHT 800

#define DO_MAP 1
#define DO_RAYS 1
#define DO_STATS 1

/* all numbers below in cm */
#define XMAX (50e-6 * 100)
#define XMIN (-50e-6 * 100)
#define YMAX (50e-6 * 100)
#define YMIN (-50e-6 * 100)
#define ZMAX (50e-6 * 100)
#define ZMIN (-50e-6 * 100)
#define STEP_SIZE 0.1

/* define simulation window ,centered on central primitive element */
#define SIM_XWIDTH (22e-6 * 100.0)
#define SIM_YWIDTH (22e-6 * 100.0)
#define SIM_ZWIDTH (22e-6 * 100.0)

/* Coordinate system ticks */
#define TICKS_MINOR (5e-6 * 100.0)
#define TICKS_MAJOR (10e-6 * 100.0)

/* model parameters */
#define RANGE_B10_LI 9.0e-6 // meter
#define RANGE_B10_HE 6.0e-6 // meter

#define RANGE_B11_HE1 3.0e-5 // meter
#define RANGE_B11_HE2 3.0e-5 // meter
#define RANGE_B11_HE3 3.0e-5 // meter

/* for this code, keep prng global */
static struct osh_rng g_rng;

static inline double urand(void) {
    return osh_rng_double(&g_rng);
}

static inline void osh_rng_unit_sphere(double out[3]) {
    /* Marsaglia (1972): sample uniformly on sphere using two uniforms */
    double u;
    double v;
    double s;

    do {
        u = 2.0 * osh_rng_double(&g_rng) - 1.0;
        v = 2.0 * osh_rng_double(&g_rng) - 1.0;
        s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);

    double const k = 2.0 * sqrt(1.0 - s);
    out[0] = u * k;
    out[1] = v * k;
    out[2] = 1.0 - 2.0 * s;
}

/* define struct which is an array of rays */
struct ray_array {
    struct ray *rays;
    double *dist; /* array of maxium distances a ray can travel */
    size_t size;  /* array length */
    char type;    /* type of reaction */
};

SDL_Color colormap[6] = {
    {255, 0, 0, 255},     // Red
    {0, 255, 0, 255},     // Green
    {100, 100, 255, 255}, // Blue
    {255, 100, 100, 255}, // Yellow
    {255, 0, 255, 255},   // Magenta
    {0, 255, 255, 255}    // Cyan
};

int random_ray(struct ray *r);
int zero_ray(struct ray *r);
int plot(struct gemca_workspace *g, int nrays, int zone);
void setRendererColor(SDL_Renderer *renderer, int zid);
static size_t material_debug_index(char const *material_name);
void drawDot(SDL_Renderer *renderer, int centerX, int centerY, int radius);
int ray_cast_statistics(struct gemca_workspace *g, int nstat);
static void draw_ray_path(SDL_Renderer *s, struct gemca_workspace *g, struct ray ray0, double max_range_cm);
static void draw_map(SDL_Renderer *s, struct gemca_workspace *g, int ndots);

int main(int argc, char *argv[]) {

    struct gemca_workspace g = {0};
    int zone = 1;

    /* Setup logger, select OSH_LOG_DEBUG for more information. */
    osh_log_init(OSH_LOG_INFO, OSH_LOG_F_NONE);

    /* Handle --version flag */
    if (argc > 1 && (strcmp(argv[1], "--version") == 0 || strcmp(argv[1], "-v") == 0)) {
        printf("bnct_sdl version %s\n", OSH_VERSION);
        return 0;
    }

    /* Handle --help flag or no arguments */
    if (argc == 1 || (argc > 1 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0))) {
        printf("Usage: %s [OPTIONS] geo.dat\n", argv[0]);
        printf("BNCT simulation and visualization tool for OpenShieldHIT\n\n");
        printf("OPTIONS:\n");
        printf("  --version, -v     Print version information\n");
        printf("  --help, -h        Show this help message\n");
        return 0;
    }

    if (argc > 1) {
        // get filename from command line and optional zone number folling the filename
        printf("Loading geometry from file: %s\n", argv[1]);
        if (argc > 2) {
            zone = atoi(argv[2]);
            printf("Using zone number: %d\n", zone);
        } else {
            zone = 1; // default to -1 if no zone number is provided
            printf("No zone number provided, defaulting to: %d\n", zone);
        }
    }

    osh_rng_init(&g_rng, OSH_RNG_TYPE_XOSHIRO256SS, 1234, 0);

    printf("----------------ffff---------------------\n");
    printf("PHASE 1: parse %s\n", argv[1]);
    fflush(stdout);
    if (osh_gemca_load(argv[1], &g) != OSH_OK) {
        fprintf(stderr, "bnct_sdl: osh_gemca_load() failed for '%s'\n", argv[1]);
        return EXIT_FAILURE;
    }

    if (g.bodies == NULL || g.zones == NULL || g.nbodies == 0 || g.nzones == 0) {
        fprintf(stderr, "bnct_sdl: geometry workspace is incomplete after osh_gemca_load()\n");
        return EXIT_FAILURE;
    }

    printf("PHASE 2: geometry loaded successfully (%llu bodies, %llu zones)\n",
           (unsigned long long) g.nbodies,
           (unsigned long long) g.nzones);
    fflush(stdout);

    osh_gemca_print_gemca(&g);

    if (plot(&g, 4000, zone) != 0) {
        fprintf(stderr, "bnct_sdl: plot() failed\n");
        return EXIT_FAILURE;
    }

    return 0;
}

/* generate a random ray in X-Z plane */
int random_ray(struct ray *r) {
    double angle;

    r->p[0] = urand() * (SIM_XWIDTH) -SIM_XWIDTH * 0.5;
    r->p[1] = 0.0;
    r->p[2] = urand() * (SIM_ZWIDTH) -SIM_ZWIDTH * 0.5;

    angle = urand() * 2.0 * M_PI;

    r->cp[0] = cos(angle);
    r->cp[1] = 0.0;
    r->cp[2] = sin(angle);

    r->system = OSH_COORD_UNIVERSE;

    return 0;
}

/* generate a random ray in X-Z plane */
int random_ray_3d(struct ray *r) {
    double cp[3]; /* direction vector */
    r->p[0] = urand() * (SIM_XWIDTH) -SIM_XWIDTH * 0.5;
    r->p[1] = urand() * (SIM_YWIDTH) -SIM_YWIDTH * 0.5;
    r->p[2] = urand() * (SIM_ZWIDTH) -SIM_ZWIDTH * 0.5;

    osh_rng_unit_sphere(cp);

    /* assign direction to first ray */
    osh_vect_copy(cp, r->cp);

    r->system = OSH_COORD_UNIVERSE;

    return 0;
}

/* generate a random ray in X-Z plane */
int random_pos(struct ray *r) {
    r->p[0] = urand() * (XMAX - XMIN) + XMIN;
    r->p[1] = 0.0;
    r->p[2] = urand() * (ZMAX - ZMIN) + ZMIN;

    r->cp[0] = 0.0;
    r->cp[1] = 0.0;
    r->cp[2] = 1.0;

    r->system = OSH_COORD_UNIVERSE;

    return 0;
}

/* generate a ray at zero along z */
int zero_ray(struct ray *r) {
    r->p[0] = 0.0;
    r->p[1] = 0.0;
    r->p[2] = 0.0;
    r->cp[0] = 0.0;
    r->cp[1] = 0.0;
    r->cp[2] = 1.0;
    r->system = OSH_COORD_UNIVERSE;
    return 0;
}

/* return two rays based on B10 reaction */
/* ray_array must be fully allocated beforehand */
int rays_b10(struct ray_array *ra) {

    double cp[3]; /* direction vector */
    size_t i;

    ra->size = 2;
    ra->dist[0] = RANGE_B10_LI;
    ra->dist[1] = RANGE_B10_HE;

    ra->rays[0].p[0] = urand() * (XMAX - XMIN) + XMIN;
    ra->rays[0].p[1] = urand() * (YMAX - YMIN) + YMIN;
    ra->rays[0].p[2] = urand() * (ZMAX - ZMIN) + ZMIN;

    for (i = 0; i < ra->size; i++) {
        ra->rays[i].system = OSH_COORD_UNIVERSE;
    }

    /* Copy the starting position to the second ray */
    osh_vect_copy(ra->rays[0].p, ra->rays[1].p);

    /* generate a random direction */
    osh_rng_unit_sphere(cp);

    /* assign direction to first ray */
    osh_vect_copy(cp, ra->rays[0].cp);

    /* assign the reversed direction to the second ray */
    osh_vect_reverse(cp, ra->rays[1].cp);

    return 1;
}

/* return three rays based on B11 reaction */
/* ray_array must be fully allocated beforehand */
int rays_b11(struct ray_array *ra) {

    double cp[3]; /* direction vector */
    size_t i;

    ra->size = 3;

    for (i = 0; i < ra->size; i++) {
        ra->rays[i].system = OSH_COORD_UNIVERSE;
    }

    ra->dist[0] = RANGE_B11_HE1;
    ra->dist[1] = RANGE_B11_HE2;
    ra->dist[2] = RANGE_B11_HE3;

    ra->rays[0].p[0] = urand() * (XMAX - XMIN) + XMIN;
    ra->rays[0].p[1] = urand() * (YMAX - YMIN) + YMIN;
    ra->rays[0].p[2] = urand() * (ZMAX - ZMIN) + ZMIN;

    /* Copy the starting position to the second ray */
    osh_vect_copy(ra->rays[0].p, ra->rays[1].p);
    osh_vect_copy(ra->rays[0].p, ra->rays[2].p);

    /* generate a random direction */
    osh_rng_unit_sphere(cp);

    /* assign direction to first ray */
    osh_vect_copy(cp, ra->rays[0].cp);

    /* generate another random direction TODO   */
    /* this is not really correct, as this particle has momentum before disintegration */
    osh_rng_unit_sphere(cp);
    osh_vect_copy(cp, ra->rays[1].cp);
    osh_vect_reverse(cp, ra->rays[2].cp);
    return 1;
}

int scanning_ray(struct ray *r) {

    static double x = XMIN;
    static double z = ZMIN;

    x += STEP_SIZE;
    if (x > XMAX) {
        x = XMIN;
        z += STEP_SIZE;
        if (z > ZMAX) {
            z = ZMIN;
        }
    }

    r->p[0] = x;
    r->p[1] = 0.0;
    r->p[2] = z;
    r->cp[0] = 0.0;
    r->cp[1] = 0.0;
    r->cp[2] = 1.0;
    r->system = OSH_COORD_UNIVERSE;
    return 0;
}

int coord2pixel(double x, double z, int *px, int *pz) {
    *px = (int) ((x - XMIN) / (XMAX - XMIN) * WINDOW_WIDTH);
    *pz = (int) ((z - ZMIN) / (ZMAX - ZMIN) * WINDOW_HEIGHT);
    return 0;
}

/* calculate a line from a given ray and a given distance */
int ray2line(struct ray *r, double d, float *x, float *z) {
    *x = r->p[0] + r->cp[0] * d;
    *z = r->p[2] + r->cp[2] * d;
    return 0;
}

int plot(struct gemca_workspace *g, int nrays, int zone) {
    int const windowHeight = WINDOW_HEIGHT;
    int const windowWidth = WINDOW_WIDTH;
    int quit = 0;

    int i;
    double x;
    double z;
    double tick_pos;
    struct ray r;

    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        fprintf(stderr, "SDL_Init(SDL_INIT_VIDEO) failed: %s\n", SDL_GetError());
        return 1;
    }

    SDL_Window *window = SDL_CreateWindow(
        "Geometry Test", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, windowWidth, windowHeight, SDL_WINDOW_SHOWN);

    if (!window) {
        fprintf(stderr, "SDL_CreateWindow() failed: %s\n", SDL_GetError());
        SDL_Quit();
        return 2;
    }

    SDL_Renderer *s = SDL_CreateRenderer(window, 0, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!s) {
        fprintf(stderr,
                "SDL_CreateRenderer(accelerated|vsync) failed: %s\nFalling back to software renderer.\n",
                SDL_GetError());
        s = SDL_CreateRenderer(window, -1, SDL_RENDERER_SOFTWARE);
    }
    if (!s) {
        fprintf(stderr, "SDL_CreateRenderer(software) failed: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 3;
    }

    SDL_Event event;

    SDL_SetRenderDrawColor(s, 0, 0, 0, 255);
    SDL_RenderClear(s);
    SDL_RenderPresent(s);

#if DO_MAP
    draw_map(s, g, 100000);
#endif

    /* Add a coordinate system */
    SDL_SetRenderDrawColor(s, 20, 50, 20, 255);
    for (i = XMIN; i < XMAX; i += 1) {
        x = (int) ((i - XMIN) / (XMAX - XMIN) * windowWidth);
        SDL_RenderDrawLine(s, x, 0, x, windowHeight);
    }
    for (i = ZMIN; i < ZMAX; i += 1) {
        z = (int) ((i - ZMIN) / (ZMAX - ZMIN) * windowHeight);
        SDL_RenderDrawLine(s, 0, z, windowWidth, z);
    }
    SDL_SetRenderDrawColor(s, 255, 255, 255, 255);
    SDL_RenderDrawLine(s, 0, windowHeight / 2, windowWidth, windowHeight / 2);
    SDL_RenderDrawLine(s, windowWidth / 2, 0, windowWidth / 2, windowHeight);

    /* Add minor and major ticks to the coordinate system */

    /* Draw Ticks */
    for (x = XMIN; x <= XMAX; x += TICKS_MAJOR) {
        tick_pos = (x - XMIN) / (XMAX - XMIN) * windowWidth;
        SDL_RenderDrawLine(s, (int) tick_pos, (windowHeight / 2) - 10, (int) tick_pos, (windowHeight / 2) + 10);
    }
    for (x = XMIN; x <= XMAX; x += TICKS_MINOR) {
        tick_pos = (x - XMIN) / (XMAX - XMIN) * windowWidth;
        SDL_RenderDrawLine(s, (int) tick_pos, (windowHeight / 2) - 5, (int) tick_pos, (windowHeight / 2) + 5);
    }
    for (z = ZMIN; z <= ZMAX; z += TICKS_MAJOR) {
        tick_pos = (z - ZMIN) / (ZMAX - ZMIN) * windowHeight;
        SDL_RenderDrawLine(s, (windowWidth / 2) - 10, (int) tick_pos, (windowWidth / 2) + 10, (int) tick_pos);
    }
    for (z = ZMIN; z <= ZMAX; z += TICKS_MINOR) {
        tick_pos = (z - ZMIN) / (ZMAX - ZMIN) * windowHeight;
        SDL_RenderDrawLine(s, (windowWidth / 2) - 5, (int) tick_pos, (windowWidth / 2) + 5, (int) tick_pos);
    }

    /* next plot a series of rays which travel trough the zones */
#if DO_RAYS
    const int batch_size = 50;
    int const max_rays = nrays;
    int total_rays = 0;

    while (total_rays < max_rays && !quit) {

        /* handle window events */
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                quit = 1;
            }
        }

        for (int k = 0; k < batch_size && total_rays < max_rays; k++) {

            random_ray(&r);

            double const range_cm = RANGE_B10_HE * 100.0;

            size_t zi0 = osh_gemca_get_zone_index(g, &r);
            if (zi0 == OSH_GEMCA_ZONE_INDEX_INVALID) {
                continue;
            }
            size_t medium0 = material_debug_index(g->zones[zi0]->material_name);

            if (medium0 != (size_t) zone) {
                continue; /* skip if we are not in a medium of interest */
            }

            int px;
            int pz;
            coord2pixel(r.p[0], r.p[2], &px, &pz);

            setRendererColor(s, (int) medium0);
            drawDot(s, px, pz, 3);

            draw_ray_path(s, g, r, range_cm);

            total_rays++;
        }
        SDL_Delay(500);

        SDL_RenderPresent(s);

        /* small delay so UI stays smooth */
        SDL_Delay(10);
    }
#endif

    // printf("Ray statistics\n");
    // i = 10;
    // while (i < 10000000) {
    //     ray_cast_statistics(g, i);
    //     i *= 10;
    // }

    SDL_RenderPresent(s);

    /* keep window up until user closes it */
    while (!quit) {
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                quit = 1;
            }
        }
        SDL_Delay(16);
    }

    SDL_DestroyRenderer(s);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}

void setRendererColor(SDL_Renderer *renderer, int zid) {
    SDL_Color color = colormap[zid % 6]; // Use zid to index into the colormap, modulo to cycle through
    SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, color.a);
}

static size_t material_debug_index(char const *material_name) {
    unsigned long value;
    char *endptr;

    if (!material_name) {
        return 0u;
    }
    if (strcmp(material_name, "blackhole") == 0 || strcmp(material_name, "vacuum") == 0) {
        return 0u;
    }

    value = strtoul(material_name, &endptr, 10);
    if (endptr != material_name && *endptr == '\0') {
        return (size_t) value;
    }

    return 0u;
}

void drawDot(SDL_Renderer *renderer, int centerX, int centerY, int radius) {
    for (int w = 0; w < radius * 2; w++) {
        for (int h = 0; h < radius * 2; h++) {
            int dx = radius - w; // Horizontal offset
            int dy = radius - h; // Vertical offset
            if ((dx * dx + dy * dy) <= (radius * radius)) {
                SDL_RenderDrawPoint(renderer, centerX + dx, centerY + dy);
            }
        }
    }
}

int ray_cast_statistics(struct gemca_workspace *g, int nstat) {
    struct ray *r;
    struct ray_array *ra;

    size_t zi;
    size_t medium;

    // char medium_name[4][16] = {"BlackHole", "InterCell", "IntraCell", "IntraNucleus"};

    double dist_rest;
    double dist;
    double hist[5];          /* per medium id */
    double hist_m[5];        /* per medium id */
    double dist_total = 0.0; /* counter for total tracklength*/
    double stat = nstat;
    double range_debug = 9e-6 * 100.0;

    int i;

    // double boron_conc[4] = {0,0,0,1};  /* relative boron conentrations */
    // double boron_conc[4] = {1,1,1,1};  /* relative boron conentrations */
    double boron_conc[4] = {0, 1, 0, 0}; /* relative boron conentrations */
    double boron_norm[4];                /* normalized boron concentrations */
    double boron_cumm[4];                /* cummulative boron concentrations */

    ra = (struct ray_array *) calloc(1, sizeof(struct ray_array));
    ra->size = 3;
    ra->type = 0;
    ra->rays = (struct ray *) calloc(ra->size, sizeof(struct ray));
    ra->dist = (double *) calloc(ra->size, sizeof(double));

    for (i = 0; i < 3; i++) {
        hist[i] = 0.0;
        hist_m[i] = 0.0;
    }

    /* normalize boron concentrations */
    for (i = 0; i < 4; i++) {
        boron_norm[i] = boron_conc[i] / (boron_conc[0] + boron_conc[1] + boron_conc[2] + boron_conc[3]);
    }

    printf("Boron concentrations: %f %f %f %f\n", boron_norm[0], boron_norm[1], boron_norm[2], boron_norm[3]);

    boron_cumm[0] = boron_norm[0];
    for (i = 1; i < 4; i++) {
        boron_cumm[i] = boron_cumm[i - 1] + boron_norm[i];
    }

    while (nstat--) {
        /* populate ray array (to be replaced  later by physics ) */
        for (size_t i = 0; i < ra->size; i++) {
            r = &ra->rays[i];
            ra->dist[i] = range_debug; /* temporary assumption on max length */
            while (1) {
                random_ray_3d(r);
                zi = osh_gemca_get_zone_index(g, r);
                if (zi == OSH_GEMCA_ZONE_INDEX_INVALID) {
                    continue;
                }
                medium = material_debug_index(g->zones[zi]->material_name);
                if (medium == 1 || medium == 2 || medium == 3) {
                    // if (medium == 3 )
                    break;
                }
                // if (urand() < boron_cumm[medium])
                //    break;
            }

            hist_m[medium]++;
        }

        /* loop over each ray */
        for (size_t i = 0; i < ra->size; i++) {
            // printf("Ray %li\n", i);
            r = &ra->rays[i];
            dist_rest = ra->dist[i];
            // printf("NEW-PARTILCE ----------------------------\n");
            // osh_transport_print_ray(r);

            while (dist_rest > 0) {
                /* get what zone and medium we are in */
                // printf("in while loop\n");
                zi = osh_gemca_get_zone_index(g, r);
                if (zi == OSH_GEMCA_ZONE_INDEX_INVALID) {
                    break;
                }
                // printf("get medium\n");
                medium = material_debug_index(g->zones[zi]->material_name);

                // printf("Zone: %li, zone_name, %s, medium: %li\n", zi, g->zones[zi]->name, medium);

                if (medium > 3) {
                    printf("Error: medium %li\n", medium);
                    exit(1);
                }

                /* get distance to next zone boundary */
                // printf("get distance\n");
                // printf("Zone: %li, medium: %li\n", zi, medium);
                // osh_transport_print_ray(r);
                dist = osh_gemca_dist(g->zones[zi], r);
                // printf("dist: %f\n", dist);

                if (dist > dist_rest) { /* stop particles travleling too far */
                    // printf("stopping particle with distance %.6e and dist_rest %.6e\n", dist, dist_rest);
                    dist = dist_rest;
                    dist_rest = 0;
                } else {
                    // printf("move particle with distance %.6e and dist_rest %.6e\n", dist, dist_rest);
                    osh_transport_move_ray(r, dist);
                    dist_rest -= dist;
                }
                hist[medium] += dist;
                dist_total += dist;
                // printf("dist: %.6e, dist_rest: %.6e\n", dist, dist_rest);

                /* move ray, and find distance to next zone */

            } /* end while ray still had distance left */
        } /* end of loop over rays */
    } /* end of loop over statistics */

    // exit(0);

    /* Normalize historgram to total distances */
    printf("stat: %e ", stat);
    for (i = 0; i < 4; i++) {
        hist[i] /= dist_total;
        /* print in 4 columns */
        printf(" %f", hist[i]);
    }
    printf("  Total distance: %f cm   check: %f sum: %f\n",
           dist_total,
           stat * range_debug * 100.0 * 3.0,
           hist[0] + hist[1] + hist[2] + hist[3]);

    printf("Medium statistics\n");
    for (i = 0; i < 4; i++) {
        printf("Medium %i: %f\n", i, hist_m[i] / stat);
    }

    free(ra->rays);
    free(ra);
    return 0;
}

static void draw_ray_path(SDL_Renderer *s, struct gemca_workspace *g, struct ray ray0, double max_range_cm) {
    struct ray r = ray0; /* work on a copy */
    double remaining = max_range_cm;
    SDL_Event dummy;
    (void) dummy;

    while (remaining > 0.0) {
        size_t zi = osh_gemca_get_zone_index(g, &r);
        if (zi == OSH_GEMCA_ZONE_INDEX_INVALID) {
            break;
        }
        size_t medium = material_debug_index(g->zones[zi]->material_name);

        double d_to_bnd = osh_gemca_dist(g->zones[zi], &r);

        /* guard: if geometry gives nonsense, avoid infinite loops */
        if (!(d_to_bnd > 0.0) || !isfinite(d_to_bnd)) {
            break;
        }

        double step = (d_to_bnd < remaining) ? d_to_bnd : remaining;

        /* compute endpoints in world coords */
        double x0 = r.p[0];
        double z0 = r.p[2];
        struct ray r2 = r;
        osh_transport_move_ray(&r2, step);
        double x1 = r2.p[0];
        double z1 = r2.p[2];

        /* map to pixels and draw */
        int px0;
        int pz0;
        int px1;
        int pz1;
        coord2pixel(x0, z0, &px0, &pz0);
        coord2pixel(x1, z1, &px1, &pz1);

        setRendererColor(s, (int) medium);
        SDL_RenderDrawLine(s, px0, pz0, px1, pz1);

        /* advance */
        r = r2;
        remaining -= step;

        /* if we exactly hit boundary, nudge forward a hair so we enter next zone
           (prevents being stuck on the boundary due to fp noise) */
        if (step == d_to_bnd && remaining > 0.0) {
            double const eps = 1e-12; /* cm; tune if needed */
            osh_transport_move_ray(&r, eps);
            remaining -= eps;
        }
    }
}

static void draw_map(SDL_Renderer *s, struct gemca_workspace *g, int ndots) {
    struct ray r;
    int i = 0;

    /* Render progressively so you see it update */
    int const present_every = 2000;

    for (i = 0; i < ndots; i++) {
        random_pos(&r);

        size_t zidx = osh_gemca_get_zone_index(g, &r);
        if (zidx == OSH_GEMCA_ZONE_INDEX_INVALID) {
            continue;
        }
        size_t medium = material_debug_index(g->zones[zidx]->material_name);

        int px;
        int pz;
        coord2pixel(r.p[0], r.p[2], &px, &pz);

        setRendererColor(s, (int) medium);
        SDL_RenderDrawPoint(s, px, pz);

        if ((i % present_every) == 0) {
            SDL_RenderPresent(s);

            /* Don’t printf every point; if you want progress: */
            // fprintf(stderr, "map: %d/%d\n", i, ndots);
        }
    }

    /* final present */
    SDL_RenderPresent(s);
}
