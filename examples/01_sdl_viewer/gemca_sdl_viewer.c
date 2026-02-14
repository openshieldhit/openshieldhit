/* simple temporary test file for playing with the new scorer */
// gcc SDLtest.c -o dsl -lSDL2 -Wall
#include <SDL.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "common/osh_const.h"
#include "common/osh_coord.h"
#include "common/osh_logger.h"
#include "gemca/osh_gemca2.h"
#include "transport/osh_transport.h"

#define WINDOW_WIDTH 800
#define WINDOW_HEIGHT 800

/* all numbers below in cm */
#define XMAX 10.0
#define XMIN -10.0
#define ZMAX 10.0
#define ZMIN -10.0
#define STEP_SIZE 0.1

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
void setRendererColor(SDL_Renderer *renderer, int zid);
void drawDot(SDL_Renderer *renderer, int centerX, int centerY, int radius);

int main(int argc, char *argv[]) {

    struct gemca_workspace g;

    /* Handle --version flag */
    if (argc > 1 && (strcmp(argv[1], "--version") == 0 || strcmp(argv[1], "-v") == 0)) {
        printf("gemca_sdl_viewer version %s\n", OSH_VERSION);
        return 0;
    }

    /* Handle --help flag or no arguments */
    if (argc == 1 || (argc > 1 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0))) {
        printf("Usage: %s [OPTIONS] geo.dat\n", argv[0]);
        printf("Geometry visualization tool for OpenShieldHIT\n\n");
        printf("OPTIONS:\n");
        printf("  --version, -v     Print version information\n");
        printf("  --help, -h        Show this help message\n");
        return 0;
    }

    printf("----------------ffff---------------------\n");
    printf("PHASE 1: parse %s\n", argv[1]);
    osh_gemca_load(argv[1], &g);
    osh_gemca_print_gemca(&g);

    plot(&g, 4);

    return 0;
}

/* generate a random ray in X-Z plane */
int random_ray(struct ray *r) {
    double angle;

    r->p[0] = ((float) rand() / RAND_MAX) * (XMAX - XMIN) + XMIN;
    r->p[1] = 0.0;
    r->p[2] = ((float) rand() / RAND_MAX) * (ZMAX - ZMIN) + ZMIN;

    angle = ((double) rand() / RAND_MAX) * 2.0 * OSH_M_PI;

    r->cp[0] = cos(angle);
    r->cp[1] = 0.0;
    r->cp[2] = sin(angle);

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
    float x, z;
    double dist;
    size_t zid; /* zone id */

    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        return 1;
    }
    SDL_Window *window = SDL_CreateWindow(
        "Geometry Test", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, windowWidth, windowHeight, SDL_WINDOW_SHOWN);
    if (window == NULL) {
        SDL_Quit();
        return 2;
    }
    // We create a renderer with hardware acceleration, we also present according with the vertical sync refresh.
    SDL_Renderer *s = SDL_CreateRenderer(window, 0, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);

    SDL_Event event;
    SDL_RenderClear(s);

    while ((!quit) && (i < ndots)) {
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                quit = 1;
            }
        }

        i++;

        // random_ray(&r);
        zero_ray(&r);
        scanning_ray(&r);
        random_ray(&r);

        zid = osh_gemca_zone(*g, r);

        /* relative position in window */
        coord2pixel(r.p[0], r.p[2], &pointLocationx, &pointLocationy);

        printf("--------------gggg----------------------(%.3f,%.3f)----- x: %i, z: %i  col:  %li - \'%s\'\n",
               r.p[0],
               r.p[2],
               pointLocationx,
               pointLocationy,
               zid,
               g->zones[zid - 1]->name);
        // Set our color for the draw functions
        setRendererColor(s, zid);
        SDL_RenderDrawPoint(s, pointLocationx, pointLocationy);
    }

    /* Add 1 cm grid */
    SDL_SetRenderDrawColor(s, 20, 50, 20, 255);
    for (i = XMIN; i < XMAX; i += 1) {
        x = (int) ((i - XMIN) / (XMAX - XMIN) * windowWidth);
        SDL_RenderDrawLine(s, x, 0, x, windowHeight);
    }
    for (i = ZMIN; i < ZMAX; i += 1) {
        z = (int) ((i - ZMIN) / (ZMAX - ZMIN) * windowHeight);
        SDL_RenderDrawLine(s, 0, z, windowWidth, z);
    }

    /* Add a coordinate system */
    SDL_SetRenderDrawColor(s, 255, 255, 255, 255);
    SDL_RenderDrawLine(s, 0, windowHeight / 2, windowWidth, windowHeight / 2);
    SDL_RenderDrawLine(s, windowWidth / 2, 0, windowWidth / 2, windowHeight);

    /* Add 1 cm scale ticks to coordinate system */
    for (i = XMIN; i < XMAX; i += 1) {
        x = (int) ((i - XMIN) / (XMAX - XMIN) * windowWidth);
        if (i % 5 == 0)
            SDL_RenderDrawLine(s, x, windowHeight / 2 - 10, x, windowHeight / 2 + 10);
        else
            SDL_RenderDrawLine(s, x, windowHeight / 2 - 5, x, windowHeight / 2 + 5);
    }
    for (i = ZMIN; i < ZMAX; i += 1) {
        z = (int) ((i - ZMIN) / (ZMAX - ZMIN) * windowHeight);
        if (i % 5 == 0)
            SDL_RenderDrawLine(s, windowWidth / 2 - 10, z, windowWidth / 2 + 10, z);
        else
            SDL_RenderDrawLine(s, windowWidth / 2 - 5, z, windowWidth / 2 + 5, z);
    }

    /* next plot a series of rays which travel trough the zones */
    for (i = 0; i < 500; i++) {
        random_ray(&r);
        // zero_ray(&r);

        /* print the ray */
        printf("Ray: (%.3f,%.3f,%.3f) (%.3f,%.3f,%.3f)\n", r.p[0], r.p[1], r.p[2], r.cp[0], r.cp[1], r.cp[2]);
        // exit(0);

        zi = osh_gemca_zone_index(*g, r);
        dist = osh_gemca_dist(g->zones[zi], &r);
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

    SDL_RenderPresent(s);

    sleep(1130);

    SDL_DestroyWindow(window);
    // We have to destroy the renderer, same as with the window.
    SDL_DestroyRenderer(s);
    SDL_Quit();

    return 1;
}

void setRendererColor(SDL_Renderer *renderer, int zid) {
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