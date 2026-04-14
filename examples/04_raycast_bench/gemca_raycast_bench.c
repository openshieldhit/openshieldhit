/* gemca_raycast_bench.c — AST vs runtime ray-casting comparison viewer
 *
 * Both panels show the same XZ plane with the same random (x, z) origins.
 * The left panel colours each point with the zone returned by the old
 * pointer-based AST lookup; the right panel uses the flat RPN runtime.
 * If both implementations agree the panels are pixel-identical; any
 * disagreement shows up immediately as a colour difference.
 *
 * After rendering the viewer also runs a standalone timing benchmark:
 * BENCH_CHUNK-sized batches of random 3-D rays are timed through the AST
 * per-ray loop and the batch runtime API.  Ray generation is excluded from
 * both timings so the numbers reflect pure geometry-query cost.
 *
 * Usage: gemca_raycast_bench [OPTIONS] geo.dat
 *   --width W      Half-extent of each view panel in cm (default 10.0)
 *   --nrays N      Benchmark ray count (default 10000000)
 *   --nrender N    Points rendered per panel (default 500000)
 *   --seed S       RNG seed (default 42)
 *   --no-display   Benchmark only, skip the SDL window
 *   --version, -v  Print version and exit
 *   --help, -h     Print this help and exit
 *
 * Press Q or close the window to quit.
 */

#include <SDL.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "common/osh_const.h"
#include "common/osh_coord.h"
#include "common/osh_logger.h"
#include "common/osh_version.h"
#include "gemca/osh_gemca2.h"
#include "gemca/runtime/osh_gemca_runtime.h"
#include "random/osh_rng.h"

/* ---- Window layout --------------------------------------------------------- */

#define PANEL_PIXELS 600
#define SEPARATOR_W 4
#define WINDOW_W (PANEL_PIXELS * 2 + SEPARATOR_W)
#define WINDOW_H PANEL_PIXELS
#define RIGHT_X (PANEL_PIXELS + SEPARATOR_W)

/* ---- Benchmark chunk size -------------------------------------------------- */

/*
 * Rays per chunk.  6 × 65536 × 8 bytes ≈ 3 MB fits comfortably in L2 cache
 * so the measurement is dominated by the zone-lookup logic, not cache misses
 * from oversized ray pools.
 */
#define BENCH_CHUNK 65536

/* ---- Render chunk size ----------------------------------------------------- */

/*
 * Zone indices accumulated before plotting — balances SDL call overhead
 * against memory used for the temporary zone arrays.
 */
// #define RENDER_CHUNK 4096
#define RENDER_CHUNK 65336 /* same as BENCH_CHUNK — no memory issue, max batch size for runtime API */

/* ---- Colormap (12 distinct colours, cycled by zone index) ------------------ */

#define N_COLORS 12
static SDL_Color colormap[N_COLORS] = {
    {220, 50, 50, 255},   /* red     */
    {50, 180, 50, 255},   /* green   */
    {60, 100, 220, 255},  /* blue    */
    {220, 180, 40, 255},  /* yellow  */
    {200, 60, 200, 255},  /* magenta */
    {40, 200, 200, 255},  /* cyan    */
    {230, 140, 40, 255},  /* orange  */
    {140, 60, 200, 255},  /* purple  */
    {180, 220, 80, 255},  /* lime    */
    {60, 180, 180, 255},  /* teal    */
    {220, 100, 140, 255}, /* pink    */
    {160, 200, 220, 255}, /* sky     */
};

/* ---- CLI options ----------------------------------------------------------- */

struct options {
    char const *geo_file;
    double half_width;
    long nrays;
    long nrender;
    uint64_t seed;
    int no_display;
};

static void print_usage(char const *prog) {
    printf("Usage: %s [OPTIONS] geo.dat\n\n", prog);
    printf("AST vs runtime zone-lookup comparison viewer and benchmark.\n\n");
    printf("  Left panel : XZ plane coloured by old AST zone lookup\n");
    printf("  Right panel: XZ plane coloured by flat RPN runtime lookup\n");
    printf("  Panels should be pixel-identical; any colour difference is a bug.\n\n");
    printf("OPTIONS:\n");
    printf("  --width W      Half-extent of each panel in cm  (default 10.0)\n");
    printf("  --nrays N      Benchmark ray count              (default 10000000)\n");
    printf("  --nrender N    Points rendered per panel        (default 500000)\n");
    printf("  --seed S       RNG seed                         (default 42)\n");
    printf("  --no-display   Benchmark only, skip SDL window\n");
    printf("  --version, -v  Print version and exit\n");
    printf("  --help, -h     Show this help and exit\n");
}

/* Returns 0 ok, -1 version/help printed (exit 0), >0 error. */
static int parse_args(int argc, char *argv[], struct options *opts) {
    int i;

    opts->geo_file = NULL;
    opts->half_width = 10.0;
    opts->nrays = 10000000L;
    opts->nrender = 500000L;
    opts->seed = 42u;
    opts->no_display = 0;

    for (i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--version") == 0 || strcmp(argv[i], "-v") == 0) {
            printf("gemca_raycast_bench version %s\n", OSH_VERSION);
            return -1;
        }
        if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            print_usage(argv[0]);
            return -1;
        }
        if (strcmp(argv[i], "--no-display") == 0) {
            opts->no_display = 1;
            continue;
        }
        if (strcmp(argv[i], "--width") == 0 && i + 1 < argc) {
            opts->half_width = atof(argv[++i]);
            if (opts->half_width <= 0.0) {
                fprintf(stderr, "error: --width must be positive\n");
                return 1;
            }
            continue;
        }
        if (strcmp(argv[i], "--nrays") == 0 && i + 1 < argc) {
            opts->nrays = atol(argv[++i]);
            if (opts->nrays <= 0) {
                fprintf(stderr, "error: --nrays must be positive\n");
                return 1;
            }
            continue;
        }
        if (strcmp(argv[i], "--nrender") == 0 && i + 1 < argc) {
            opts->nrender = atol(argv[++i]);
            if (opts->nrender <= 0) {
                fprintf(stderr, "error: --nrender must be positive\n");
                return 1;
            }
            continue;
        }
        if (strcmp(argv[i], "--seed") == 0 && i + 1 < argc) {
            opts->seed = (uint64_t) atol(argv[++i]);
            continue;
        }
        if (argv[i][0] != '-') {
            opts->geo_file = argv[i];
            continue;
        }
        fprintf(stderr, "error: unknown option '%s'\n", argv[i]);
        return 1;
    }

    if (!opts->geo_file) {
        print_usage(argv[0]);
        return 1;
    }
    return 0;
}

/* ---- Coordinate mapping ---------------------------------------------------- */

static inline int coord_to_pixel(double val, double half) {
    int px = (int) ((val + half) / (2.0 * half) * (double) PANEL_PIXELS);
    if (px < 0) {
        px = 0;
    }
    if (px >= PANEL_PIXELS) {
        px = PANEL_PIXELS - 1;
    }
    return px;
}

/* ---- Ray generation -------------------------------------------------------- */

/*
 * Fill SoA arrays for `chunk` rays: random 3-D positions in [-half, half]^3
 * with cheap random directions.
 *
 * Zone membership is purely positional — direction has no effect on which zone
 * a point falls in.  Instead of generating a proper isotropic direction (which
 * requires sincos and costs ~12% of profiling time), we draw two flat uniforms
 * and map them to [-1, 1] for the X and Z components; Y is fixed at zero.
 * The direction is not normalised — that is only needed for distance queries.
 *
 * Three batch RNG calls fill positions; two fill direction components.
 * The transform loop is then pure linear arithmetic, fully auto-vectorisable.
 */
static void generate_chunk(struct osh_rng *rng,
                           double half,
                           long chunk,
                           double *bx,
                           double *by,
                           double *bz,
                           double *bux,
                           double *buy,
                           double *buz) {
    long i;
    int n = (int) chunk; /* BENCH_CHUNK <= 65536, safe to cast */

    /* Batch-generate raw uniforms [0, 1) for positions and two direction components. */
    osh_rng_double_vec(rng, bx,  n);
    osh_rng_double_vec(rng, by,  n);
    osh_rng_double_vec(rng, bz,  n);
    osh_rng_double_vec(rng, bux, n); /* X direction raw → [-1, 1] */
    osh_rng_double_vec(rng, buz, n); /* Z direction raw → [-1, 1] */

    /* Transform in-place — no trig, fully vectorisable. */
    for (i = 0; i < chunk; ++i) {
        bx[i]  = bx[i]  * 2.0 * half - half;
        by[i]  = by[i]  * 2.0 * half - half;
        bz[i]  = bz[i]  * 2.0 * half - half;
        bux[i] = bux[i] * 2.0 - 1.0;
        buy[i] = 0.0;
        buz[i] = buz[i] * 2.0 - 1.0;
    }
}

/* ---- Benchmark ------------------------------------------------------------- */

/*
 * Time zone lookup for `nrays` rays through both APIs.
 *
 * Per chunk:
 *   1. generate_chunk() fills the SoA arrays (not timed).
 *   2. AST per-ray loop is timed.
 *   3. Batch runtime call is timed on the same rays.
 *
 * A volatile sink accumulates XOR'd zone indices to prevent dead-code
 * elimination by the compiler.
 */
static void run_benchmark(struct gemca_workspace *g,
                          struct gemca_runtime const *rt,
                          double *bx,
                          double *by,
                          double *bz,
                          double *bux,
                          double *buy,
                          double *buz,
                          size_t *zone_out,
                          double half,
                          long nrays,
                          uint64_t seed) {
    struct osh_rng rng;
    struct ray ray;
    long done;
    long chunk;
    long i;
    clock_t t0;
    clock_t t1;
    clock_t ticks_ast;
    clock_t ticks_batch;
    size_t volatile sink;
    double elapsed_ast;
    double elapsed_batch;

    printf("\n--- Benchmark: %ld rays, chunk %d ---\n", nrays, BENCH_CHUNK);
    printf("    Geometry: %zu zone(s), %zu body/bodies\n\n", g->nzones, g->nbodies);

    ray.system = OSH_COORD_UNIVERSE;
    ticks_ast = 0;
    ticks_batch = 0;
    sink = 0u;

    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, seed, 10u /* independent stream */);
    done = 0;
    while (done < nrays) {
        chunk = nrays - done;
        if (chunk > BENCH_CHUNK) {
            chunk = BENCH_CHUNK;
        }

        /* Step 1 — generate rays (excluded from both timings) */
        generate_chunk(&rng, half, chunk, bx, by, bz, bux, buy, buz);

        /* Step 2 — AST per-ray */
        t0 = clock();
        for (i = 0; i < chunk; ++i) {
            ray.p[0] = bx[i];
            ray.p[1] = by[i];
            ray.p[2] = bz[i];
            ray.cp[0] = bux[i];
            ray.cp[1] = buy[i];
            ray.cp[2] = buz[i];
            sink ^= osh_gemca_get_zone_index(g, &ray);
        }
        t1 = clock();
        ticks_ast += t1 - t0;

        /* Step 3 — runtime batch on the same rays */
        t0 = clock();
        osh_gemca_runtime_get_zone_batch(rt, bx, by, bz, bux, buy, buz, (size_t) chunk, zone_out);
        t1 = clock();
        ticks_batch += t1 - t0;

        for (i = 0; i < chunk; ++i) {
            sink ^= zone_out[i];
        }

        done += chunk;
    }
    (void) sink;

    elapsed_ast = (double) ticks_ast / (double) CLOCKS_PER_SEC;
    elapsed_batch = (double) ticks_batch / (double) CLOCKS_PER_SEC;

    printf("    AST (per-ray):    %8.3f s  —  %6.2f M rays/s\n", elapsed_ast, (double) nrays / elapsed_ast / 1e6);
    printf("    Runtime (batch):  %8.3f s  —  %6.2f M rays/s\n", elapsed_batch, (double) nrays / elapsed_batch / 1e6);
    if (elapsed_batch > 0.0) {
        printf("    Speedup:          %6.2fx\n\n", elapsed_ast / elapsed_batch);
    }
}

/* ---- SDL rendering --------------------------------------------------------- */

static void draw_grid(SDL_Renderer *renderer, int x_off, double half) {
    int px;
    double v;

    SDL_SetRenderDrawColor(renderer, 25, 45, 25, 255);
    for (v = -half; v <= half; v += 1.0) {
        px = coord_to_pixel(v, half);
        SDL_RenderDrawLine(renderer, x_off + px, 0, x_off + px, WINDOW_H);
        SDL_RenderDrawLine(renderer, x_off, px, x_off + PANEL_PIXELS, px);
    }

    /* Bright axis crosshairs */
    SDL_SetRenderDrawColor(renderer, 180, 180, 180, 255);
    px = coord_to_pixel(0.0, half);
    SDL_RenderDrawLine(renderer, x_off + px, 0, x_off + px, WINDOW_H);
    SDL_RenderDrawLine(renderer, x_off, px, x_off + PANEL_PIXELS, px);
}

static void draw_separator(SDL_Renderer *renderer) {
    SDL_Rect sep;
    sep.x = PANEL_PIXELS;
    sep.y = 0;
    sep.w = SEPARATOR_W;
    sep.h = WINDOW_H;
    SDL_SetRenderDrawColor(renderer, 70, 70, 70, 255);
    SDL_RenderFillRect(renderer, &sep);
}

/*
 * Render the comparison view: both panels show the same XZ (y = 0) origins.
 *
 * For each chunk of RENDER_CHUNK points:
 *   - positions are generated once (same x, z for both panels)
 *   - left panel:  colour from osh_gemca_get_zone_index (AST)
 *   - right panel: colour from osh_gemca_runtime_get_zone_batch (runtime)
 *
 * If both implementations are correct the panels are pixel-identical.
 * A zone-index mismatch shows as a different colour on the right panel.
 *
 * Direction is fixed along +Z.  Zone membership is position-only, so direction
 * has no effect on the comparison.
 */
static void render_comparison(SDL_Renderer *renderer,
                              struct gemca_workspace *g,
                              struct gemca_runtime const *rt,
                              double half,
                              struct osh_rng *rng,
                              long npoints) {
    /* Stack-allocated chunk buffers — RENDER_CHUNK is small enough. */
    double rx[RENDER_CHUNK];
    double rz[RENDER_CHUNK];
    double ry_fixed[RENDER_CHUNK]; /* y = 0 for all */
    double rux[RENDER_CHUNK];
    double ruy[RENDER_CHUNK];
    double ruz[RENDER_CHUNK];
    size_t zone_rt[RENDER_CHUNK];

    struct ray ray_ast;
    SDL_Color col;
    size_t zone_idx_ast;
    size_t zone_idx_rt;
    long done;
    long chunk;
    long i;
    int px;
    int pz;
    long mismatches = 0;

    ray_ast.system = OSH_COORD_UNIVERSE;
    ray_ast.p[1] = 0.0; /* XZ plane: y = 0 */
    ray_ast.cp[0] = 0.0;
    ray_ast.cp[1] = 0.0;
    ray_ast.cp[2] = 1.0; /* direction along +Z — irrelevant for membership */

    /* Pre-fill constant arrays */
    for (i = 0; i < RENDER_CHUNK; ++i) {
        ry_fixed[i] = 0.0;
        rux[i] = 0.0;
        ruy[i] = 0.0;
        ruz[i] = 1.0;
    }

    done = 0;
    while (done < npoints) {
        chunk = npoints - done;
        if (chunk > RENDER_CHUNK) {
            chunk = RENDER_CHUNK;
        }

        /* Generate random (x, z) positions in the XZ plane using batch RNG. */
        osh_rng_double_vec(rng, rx, (int) chunk);
        osh_rng_double_vec(rng, rz, (int) chunk);
        for (i = 0; i < chunk; ++i) {
            rx[i] = rx[i] * 2.0 * half - half;
            rz[i] = rz[i] * 2.0 * half - half;
        }

        /* Batch runtime query for the whole chunk */
        osh_gemca_runtime_get_zone_batch(rt, rx, ry_fixed, rz, rux, ruy, ruz, (size_t) chunk, zone_rt);

        /* Per-point: AST query + plot both panels */
        for (i = 0; i < chunk; ++i) {
            ray_ast.p[0] = rx[i];
            ray_ast.p[2] = rz[i];

            zone_idx_ast = osh_gemca_get_zone_index(g, &ray_ast);
            zone_idx_rt = zone_rt[i];

            px = coord_to_pixel(rx[i], half);
            pz = coord_to_pixel(rz[i], half);

            /* Left panel — AST */
            if (zone_idx_ast != OSH_GEMCA_ZONE_INDEX_INVALID) {
                col = colormap[zone_idx_ast % N_COLORS];
                SDL_SetRenderDrawColor(renderer, col.r, col.g, col.b, col.a);
                SDL_RenderDrawPoint(renderer, px, pz);
            }

            /* Right panel — runtime */
            if (zone_idx_rt != OSH_GEMCA_ZONE_INDEX_INVALID) {
                col = colormap[zone_idx_rt % N_COLORS];
                SDL_SetRenderDrawColor(renderer, col.r, col.g, col.b, col.a);
                SDL_RenderDrawPoint(renderer, RIGHT_X + px, pz);
            }

            if (zone_idx_ast != zone_idx_rt) {
                ++mismatches;
            }
        }

        done += chunk;
    }

    printf("Render complete: %ld points, %ld mismatch(es) between AST and runtime.\n", npoints, mismatches);
    if (mismatches > 0) {
        printf("WARNING: mismatches detected — check zone boundary logic.\n");
    }
}

/* ---- Main ------------------------------------------------------------------ */

int main(int argc, char *argv[]) {
    struct options opts;
    struct gemca_workspace *g = NULL;
    struct gemca_runtime rt;
    struct osh_rng rng;
    double *bx = NULL;
    double *by = NULL;
    double *bz = NULL;
    double *bux = NULL;
    double *buy = NULL;
    double *buz = NULL;
    size_t *zone_out = NULL;
    SDL_Window *window = NULL;
    SDL_Renderer *renderer = NULL;
    SDL_Event event;
    char title[320];
    int rc;
    int quit;

    memset(&rt, 0, sizeof(rt));

    osh_log_init(OSH_LOG_INFO, OSH_LOG_F_NONE);

    rc = parse_args(argc, argv, &opts);
    if (rc != 0) {
        return (rc < 0) ? 0 : 1;
    }

    /* ---- Load geometry ---------------------------------------------------- */

    if (osh_gemca_workspace_init(&g) != OSH_OK) {
        fprintf(stderr, "error: failed to allocate gemca workspace\n");
        return 1;
    }
    if (osh_gemca_load(opts.geo_file, g) != OSH_OK) {
        fprintf(stderr, "error: failed to load geometry from '%s'\n", opts.geo_file);
        osh_gemca_workspace_free(g);
        return 1;
    }
    osh_gemca_print_gemca(g);

    /* ---- Compile runtime -------------------------------------------------- */

    if (osh_gemca_runtime_setup(g, &rt) != OSH_OK) {
        fprintf(stderr, "error: failed to compile gemca runtime\n");
        osh_gemca_workspace_free(g);
        return 1;
    }

    /* ---- Benchmark scratch (heap) ----------------------------------------- */

    bx = malloc((size_t) BENCH_CHUNK * sizeof(double));
    by = malloc((size_t) BENCH_CHUNK * sizeof(double));
    bz = malloc((size_t) BENCH_CHUNK * sizeof(double));
    bux = malloc((size_t) BENCH_CHUNK * sizeof(double));
    buy = malloc((size_t) BENCH_CHUNK * sizeof(double));
    buz = malloc((size_t) BENCH_CHUNK * sizeof(double));
    zone_out = malloc((size_t) BENCH_CHUNK * sizeof(size_t));

    if (!bx || !by || !bz || !bux || !buy || !buz || !zone_out) {
        fprintf(stderr, "error: benchmark allocation failed\n");
        rc = 1;
        goto cleanup;
    }

    /* ---- Benchmark --------------------------------------------------------- */

    run_benchmark(g, &rt, bx, by, bz, bux, buy, buz, zone_out, opts.half_width, opts.nrays, opts.seed);

    if (opts.no_display) {
        goto cleanup;
    }

    /* ---- SDL visualisation ------------------------------------------------- */

    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        fprintf(stderr, "SDL_Init failed: %s\n", SDL_GetError());
        rc = 1;
        goto cleanup;
    }

    snprintf(title,
             sizeof(title),
             "Left: AST  |  Right: runtime  |  XZ plane (y=0)  |"
             "  half=%.1f cm  |  %zu zones  |  %s",
             opts.half_width,
             rt.nzones,
             opts.geo_file);

    window =
        SDL_CreateWindow(title, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WINDOW_W, WINDOW_H, SDL_WINDOW_SHOWN);
    if (!window) {
        fprintf(stderr, "SDL_CreateWindow failed: %s\n", SDL_GetError());
        rc = 1;
        goto sdl_cleanup;
    }

    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!renderer) {
        fprintf(stderr, "SDL_CreateRenderer failed: %s\n", SDL_GetError());
        rc = 1;
        goto sdl_cleanup;
    }

    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);

    /* Render stream 2 — independent from benchmark stream 10 */
    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, opts.seed, 2u);

    printf("Rendering %ld points per panel (XZ plane, y=0)...\n", opts.nrender);
    render_comparison(renderer, g, &rt, opts.half_width, &rng, opts.nrender);

    draw_grid(renderer, 0, opts.half_width);
    draw_grid(renderer, RIGHT_X, opts.half_width);
    draw_separator(renderer);

    SDL_RenderPresent(renderer);
    printf("Left=AST  Right=runtime  |  Press Q or close window to exit.\n");

    quit = 0;
    while (!quit) {
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                quit = 1;
            }
            if (event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_q) {
                quit = 1;
            }
        }
        SDL_Delay(16);
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);

sdl_cleanup:
    SDL_Quit();

cleanup:
    free(bx);
    free(by);
    free(bz);
    free(bux);
    free(buy);
    free(buz);
    free(zone_out);
    osh_gemca_runtime_free(&rt);
    osh_gemca_workspace_free(g);
    return rc;
}
