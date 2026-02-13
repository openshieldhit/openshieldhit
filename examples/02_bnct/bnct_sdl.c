/* simple temporary test file for playing with the new scorer */
// gcc SDLtest.c -o dsl -lSDL2 -Wall
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <SDL2/SDL.h>
#include <math.h>

#include "common/sh_logger.h"
#include "common/sh_coord.h"
#include "common/sh_vect.h"
#include "transport/sh_transport.h"
#include "random/sh_random.h"

#include "gemca/sh_gemca2.h"

#define WINDOW_WIDTH 800
#define WINDOW_HEIGHT 800

/* all numbers below in cm */
#define XMAX  50e-6 * 100
#define XMIN -50e-6 * 100
#define YMAX  50e-6 * 100
#define YMIN -50e-6 * 100
#define ZMAX  50e-6 * 100
#define ZMIN -50e-6 * 100
#define STEP_SIZE 0.1

/* define simulation window ,centered on central primitive element */
#define SIM_XWIDTH  22e-6 * 100.0
#define SIM_YWIDTH  22e-6 * 100.0
#define SIM_ZWIDTH  22e-6 * 100.0

/* Coordinate system ticks */
#define TICKS_MINOR  5e-6 * 100.0
#define TICKS_MAJOR  10e-6 * 100.0

/* model parameters */
#define RANGE_B10_LI 9.0e-6 // meter
#define RANGE_B10_HE 6.0e-6 // meter

#define RANGE_B11_HE1 3.0e-5 // meter
#define RANGE_B11_HE2 3.0e-5 // meter
#define RANGE_B11_HE3 3.0e-5 // meter

/* define struct which is an array of rays */
struct ray_array {
    struct ray *rays;
    double *dist;  /* array of maxium distances a ray can travel */
    size_t size; /* array length */
    char type;  /* type of reaction */
};


SDL_Color colormap[6] = {
    {255, 0, 0, 255},   // Red
    {0, 255, 0, 255},   // Green
    {0, 0, 255, 255},   // Blue
    {255, 255, 0, 255}, // Yellow
    {255, 0, 255, 255}, // Magenta
    {0, 255, 255, 255}  // Cyan
};

int random_ray(struct ray *r);
int zero_ray(struct ray *r);
int plot(struct gemca_workspace *g, int ndots);
void setRendererColor(SDL_Renderer* renderer, int zid);
void drawDot(SDL_Renderer *renderer, int centerX, int centerY, int radius);
int ray_cast_statistics(struct gemca_workspace *g, int nstat);
int test_ray_trouble(struct ray *r);

int main(int argc, char *argv[]) {

    struct gemca_workspace g;

    sh_random_init(1234);

    if (argc == 1) {
        printf("Usage: %s geo.dat\n", argv[0]);
        exit(0);
    }

    printf("----------------ffff---------------------\n");
    printf("PHASE 1: parse %s\n", argv[1]);
    sh_gemca_load(argv[1], &g);
    sh_gemca_print_gemca(&g);

    plot(&g, 40000);

    return 0;
}


/* generate a random ray in X-Z plane */
int random_ray(struct ray *r) {
    double angle;

    r->p[0] = sh_random() * (SIM_XWIDTH)  -SIM_XWIDTH * 0.5;
    r->p[1] = 0.0;
    r->p[2] = sh_random() * (SIM_ZWIDTH) -SIM_ZWIDTH * 0.5;

    angle = sh_random() * 2.0 * M_PI;

    r->cp[0] = cos(angle);
    r->cp[1] = 0.0;
    r->cp[2] = sin(angle);

    r->system = SH_COORD_UNIVERSE;

    return 0;
}


/* generate a random ray in X-Z plane */
int random_ray_3d(struct ray *r) {
    double cp[3]; /* direction vector */
    r->p[0] = sh_random() * (SIM_XWIDTH) -SIM_XWIDTH * 0.5;
    r->p[1] = sh_random() * (SIM_YWIDTH) -SIM_YWIDTH * 0.5;
    r->p[2] = sh_random() * (SIM_ZWIDTH) -SIM_ZWIDTH * 0.5;

    sh_random_sphere(cp);

    /* assign direction to first ray */
    sh_vect_copy(cp, r->cp);

    r->system = SH_COORD_UNIVERSE;

    return 0;
}


/* generate a random ray in X-Z plane */
int random_pos(struct ray *r) {
    r->p[0] = sh_random() * (XMAX - XMIN) + XMIN;
    r->p[1] = 0.0;
    r->p[2] = sh_random() * (ZMAX - ZMIN) + ZMIN;

    r->cp[0] = 0.0;
    r->cp[1] = 0.0;
    r->cp[2] = 1.0;

    r->system = SH_COORD_UNIVERSE;

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
    r->system = SH_COORD_UNIVERSE;
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

    ra->rays[0].p[0] = sh_random() * (XMAX - XMIN) + XMIN;
    ra->rays[0].p[1] = sh_random() * (YMAX - YMIN) + YMIN;
    ra->rays[0].p[2] = sh_random() * (ZMAX - ZMIN) + ZMIN;

    for(i = 0; i < ra->size; i++) {
        ra->rays[i].system = SH_COORD_UNIVERSE;
    }

    /* Copy the starting position to the second ray */
    sh_vect_copy(ra->rays[0].p, ra->rays[1].p);

    /* generate a random direction */
    sh_random_sphere(cp);

    /* assign direction to first ray */
    sh_vect_copy(cp, ra->rays[0].cp);

    /* assign the reversed direction to the second ray */
    sh_vect_reverse(cp, ra->rays[1].cp);

    return 1;
}


/* return three rays based on B11 reaction */
/* ray_array must be fully allocated beforehand */
int rays_b11(struct ray_array *ra) {

    double cp[3]; /* direction vector */
    size_t i;

    ra->size = 3;

    for(i = 0; i < ra->size; i++) {
        ra->rays[i].system = SH_COORD_UNIVERSE;
    }

    ra->dist[0] = RANGE_B11_HE1;
    ra->dist[1] = RANGE_B11_HE2;
    ra->dist[2] = RANGE_B11_HE3;


    ra->rays[0].p[0] = sh_random() * (XMAX - XMIN) + XMIN;
    ra->rays[0].p[1] = sh_random() * (YMAX - YMIN) + YMIN;
    ra->rays[0].p[2] = sh_random() * (ZMAX - ZMIN) + ZMIN;

    /* Copy the starting position to the second ray */
    sh_vect_copy(ra->rays[0].p, ra->rays[1].p);
    sh_vect_copy(ra->rays[0].p, ra->rays[2].p);

    /* generate a random direction */
    sh_random_sphere(cp);

    /* assign direction to first ray */
    sh_vect_copy(cp, ra->rays[0].cp);


    /* generate another random direction TODO   */
    /* this is not really correct, as this particle has momentum before disintegration */
    sh_random_sphere(cp);
    sh_vect_copy(cp, ra->rays[1].cp);
    sh_vect_reverse(cp, ra->rays[2].cp);
    return 1;
}


int test_ray_trouble(struct ray *r) {
    r->p[0] = 0.000467384;
    r->p[1] = 0.0;
    r->p[2] = -0.001259760;
    r->cp[0] = 0.075920345;
    r->cp[1] = 0.0;
    r->cp[2] = -0.997113886;
    r->system = SH_COORD_UNIVERSE;
    return 0;
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
    r->system = SH_COORD_UNIVERSE;
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


int plot(struct gemca_workspace *g, int ndots) {
    const int windowHeight = WINDOW_HEIGHT;
    const int windowWidth = WINDOW_WIDTH;
    int pointLocationx;
    int pointLocationy;
    int pointLocationx2;
    int pointLocationy2;
    int i = 0;
    int quit = 0;
    int zi;


    struct ray r;
    float x, z, _f;
    double dist;
    size_t zid;     /* zone id */
    size_t zidx;     /* zone id */
    size_t medium;

    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        return 1;
    }
    SDL_Window *window = SDL_CreateWindow("Geometry Test",
                                          SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, windowWidth,
                                          windowHeight, SDL_WINDOW_SHOWN);
    if (window == NULL) {
        SDL_Quit();
        return 2;
    }
    // We create a renderer with hardware acceleration, we also present according with the vertical sync refresh.
    SDL_Renderer *s = SDL_CreateRenderer(window, 0, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);

    SDL_Event event;
    SDL_RenderClear(s);

    /* make a MC map of the canvas */
    while ((!quit) && (i < ndots)) {
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                quit = 1;
            }
        }

        i++;

        // random_ray(&r);
        // zero_ray(&r);
        // scanning_ray(&r);
        random_pos(&r);

        zid = sh_gemca_zone(*g, r);
        zidx = sh_gemca_zone_index(*g, r);
        medium = g->zones[zidx]->medium;


        /* relative position in window */
        coord2pixel(r.p[0], r.p[2], &pointLocationx, &pointLocationy);

        printf("--------------gggg----------------------(%.3f,%.3f)----- x: %i, z: %i  col:  %li - \'%s\'\n",
               r.p[0], r.p[2],
               pointLocationx, pointLocationy, zid, g->zones[zid-1]->name);
        // Set our color for the draw functions
        setRendererColor(s, medium);
        SDL_RenderDrawPoint(s, pointLocationx, pointLocationy);
    }

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
        _f = (int)((x - XMIN) / (XMAX - XMIN) * windowWidth);
        SDL_RenderDrawLine(s, _f, windowHeight / 2 - 10, _f, windowHeight / 2 + 10); // Major tick
    }
    for (x = XMIN; x <= XMAX; x += TICKS_MINOR) {
        _f = (int)((x - XMIN) / (XMAX - XMIN) * windowWidth);
        SDL_RenderDrawLine(s, _f, windowHeight / 2 - 5, _f, windowHeight / 2 + 5); // Major tick
    }
    for (z = ZMIN; z <= ZMAX; z += TICKS_MAJOR) {
        _f = (int)((z - ZMIN) / (ZMAX - ZMIN) * windowHeight);
        SDL_RenderDrawLine(s, windowWidth / 2 - 10, _f, windowWidth / 2 + 10, _f); // Major tick
    }
    for (z = ZMIN; z <= ZMAX; z += TICKS_MINOR) {
        _f = (int)((z - ZMIN) / (ZMAX - ZMIN) * windowHeight);
        SDL_RenderDrawLine(s, windowWidth / 2 - 5, _f, windowWidth / 2 + 5, _f); // Major tick
    }


    /* next plot a series of rays which travel trough the zones */
    for (i = 0; i < 50; i++) {
        random_ray(&r);
        //zero_ray(&r);

        /* print the ray */
        printf("Ray: (%.3f,%.3f,%.3f) (%.3f,%.3f,%.3f)\n", r.p[0], r.p[1], r.p[2], r.cp[0], r.cp[1], r.cp[2]);
        // exit(0);

        zi = sh_gemca_zone_index(*g, r);
        dist = sh_gemca_dist(g->zones[zi], &r);
        printf("dist: %f\n", dist);
        if (dist > 20.0)
            dist = 20.0;

        ray2line(&r, dist, &x, &z);

        /* get coordinates of line */
        coord2pixel(r.p[0], r.p[2], &pointLocationx, &pointLocationy);
        coord2pixel(x, z, &pointLocationx2, &pointLocationy2);

        /* draw each line segment with a new color */
        setRendererColor(s, zi);
        SDL_RenderDrawLine(s, pointLocationx, pointLocationy, pointLocationx2, pointLocationy2);
        /* and add a thick dot for the start position of the ray */
        drawDot(s, pointLocationx, pointLocationy, 3);
    }

    printf("Ray statistics\n");
    i = 10;
    while(i < 10000000) {
        ray_cast_statistics(g, i);
        i *= 10;
    }

    exit(0);
    SDL_RenderPresent(s);

    sleep(1130);

    SDL_DestroyWindow(window);
    // We have to destroy the renderer, same as with the window.
    SDL_DestroyRenderer(s);
    SDL_Quit();

    return 1;
}


void setRendererColor(SDL_Renderer* renderer, int zid) {
    SDL_Color color = colormap[zid % 6]; // Use zid to index into the colormap, modulo to cycle through
    SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, color.a);
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
    double hist[5]; /* per medium id */
    double hist_m[5]; /* per medium id */
    double dist_total = 0.0; /* counter for total tracklength*/
    double stat = nstat;
    double range_debug = 9e-6 * 100.0;

    int i;

    // double boron_conc[4] = {0,0,0,1};  /* relative boron conentrations */
    //double boron_conc[4] = {1,1,1,1};  /* relative boron conentrations */
    double boron_conc[4] = {0,1,0,0};  /* relative boron conentrations */
    double boron_norm[4];  /* normalized boron concentrations */
    double boron_cumm[4];  /* cummulative boron concentrations */

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


    while(nstat--) {
        /* populate ray array (to be replaced  later by physics ) */
        for (size_t i = 0; i < ra->size; i++) {
            r = &ra->rays[i];
            ra->dist[i] = range_debug;  /* temporary assumption on max length */
            while(1) {
                random_ray_3d(r);
                zi = sh_gemca_zone_index(*g, *r);
                medium = g->zones[zi]->medium; /* check medium */
                if (medium == 1 || medium == 2 || medium == 3)
                    // if (medium == 3 )
                    break;
                // if (sh_random() < boron_cumm[medium])
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
            // sh_transport_print_ray(r);

            while(dist_rest > 0) {
                /* get what zone and medium we are in */
                // printf("in while loop\n");
                zi = sh_gemca_zone_index(*g, *r);
                // printf("get medium\n");
                medium = g->zones[zi]->medium;

                // printf("Zone: %li, zone_name, %s, medium: %li\n", zi, g->zones[zi]->name, medium);

                if (medium > 3) {
                    printf("Error: medium %li\n", medium);
                    exit(1);
                }

                /* get distance to next zone boundary */
                // printf("get distance\n");
                // printf("Zone: %li, medium: %li\n", zi, medium);
                // sh_transport_print_ray(r);
                dist = sh_gemca_dist(g->zones[zi], r);
                // printf("dist: %f\n", dist);

                if (dist > dist_rest) { /* stop particles travleling too far */
                    // printf("stopping particle with distance %.6e and dist_rest %.6e\n", dist, dist_rest);
                    dist = dist_rest;
                    dist_rest = 0;
                } else {
                    // printf("move particle with distance %.6e and dist_rest %.6e\n", dist, dist_rest);
                    sh_transport_move_ray(r, dist);
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
    printf("  Total distance: %f cm   check: %f sum: %f\n", dist_total, stat * range_debug * 100.0*3.0,
           hist[0] + hist[1] + hist[2] + hist[3]);

    printf("Medium statistics\n");
    for (i = 0; i < 4; i++) {
        printf("Medium %i: %f\n", i, hist_m[i] / stat);
    }

    free(ra->rays);
    free(ra);
    return 0;
}