#include <stdlib.h>
#include <string.h>

#include "common/osh_logger.h"
#include "gemca/osh_gemca2.h"
#include "gemca/osh_gemca2_calc_body.h"
#include "gemca/osh_gemca2_internal.h"
#include "gemca/prepare/osh_gemca_body_prepare.h"
#include "gemca/prepare/osh_gemca_zone_compile.h"
#include "openshieldhit/geometry.h"

/* ---- Internal helpers ---------------------------------------------------- */

static void geometry_body_free_fields(struct osh_geometry_body *b) {
    if (!b) {
        return;
    }
    free(b->name);
    free(b->a);
    b->name = NULL;
    b->a = NULL;
}

static void geometry_zone_free_fields(struct osh_geometry_zone *z) {
    if (!z) {
        return;
    }
    free(z->name);
    free(z->material_name);
    free(z->expr);
    z->name = NULL;
    z->material_name = NULL;
    z->expr = NULL;
}

/* ---- Lifecycle ----------------------------------------------------------- */

int osh_geometry_workspace_create(struct osh_geometry_workspace **ws_out) {
    struct osh_geometry_workspace *ws;

    if (!ws_out) {
        return -1;
    }
    ws = (struct osh_geometry_workspace *) calloc(1, sizeof(*ws));
    if (!ws) {
        return -1;
    }
    *ws_out = ws;
    return 0;
}

int osh_geometry_workspace_prepare(struct osh_geometry_workspace *ws) {
    size_t i;
    struct osh_gemca_prepared *gemca = NULL;

    if (!ws) {
        return -1;
    }
    if (ws->prepared) {
        osh_error("%s", "geometry: osh_geometry_workspace_prepare called on an already-prepared workspace");
        return -1;
    }

    /* ---- 1. Allocate internal gemca workspace -------------------------------- */

    gemca = (struct osh_gemca_prepared *) calloc(1, sizeof(*gemca));
    if (!gemca) {
        return -1;
    }

    /* ---- 2. Build internal body objects from cold data ----------------------- */

    gemca->nbodies = ws->nbodies;

    if (ws->nbodies > 0u) {
        gemca->bodies = (struct body **) calloc(ws->nbodies, sizeof(struct body *));
        if (!gemca->bodies) {
            goto fail;
        }
    }

    for (i = 0u; i < ws->nbodies; ++i) {
        struct osh_geometry_body const *cb = &ws->bodies[i];
        struct body *b;

        if (osh_gemca_body_init(&b) != OSH_OK) {
            goto fail;
        }
        gemca->bodies[i] = b;

        b->type = cb->type;
        b->coord = cb->coord;
        b->na = cb->na;

        if (cb->name) {
            b->name = strdup(cb->name);
            if (!b->name) {
                goto fail;
            }
        }

        if (cb->na > 0 && cb->a) {
            b->a = (double *) calloc((size_t) cb->na, sizeof(double));
            if (!b->a) {
                goto fail;
            }
            memcpy(b->a, cb->a, (size_t) cb->na * sizeof(double));
        }
    }

    /* Compute surfaces and transformation matrices for all bodies. */
    if (osh_gemca_body_setup(gemca) != OSH_OK) {
        goto fail;
    }

    /* ---- 3. Build internal zone objects and compile their ASTs --------------- */

    gemca->nzones = ws->nzones;

    if (ws->nzones > 0u) {
        gemca->zones = (struct zone **) calloc(ws->nzones, sizeof(struct zone *));
        if (!gemca->zones) {
            goto fail;
        }
    }

    for (i = 0u; i < ws->nzones; ++i) {
        struct osh_geometry_zone const *cz = &ws->zones[i];
        struct zone *z;

        if (osh_gemca_zone_init(&z) != OSH_OK) {
            goto fail;
        }
        gemca->zones[i] = z;

        if (cz->name) {
            z->name = strdup(cz->name);
            if (!z->name) {
                goto fail;
            }
        }

        if (cz->material_name) {
            z->material_name = strdup(cz->material_name);
            if (!z->material_name) {
                goto fail;
            }
        }

        /* Compile the boolean expression into a cgnode AST.  Bodies must be
         * set up already (step 2 above) because the tokeniser resolves body
         * names to pointers at this stage. */
        if (cz->expr && cz->expr[0]) {
            if (osh_gemca_zone_compile_expr(z, cz->expr, gemca) != OSH_OK) {
                goto fail;
            }
        }
    }

    ws->prepared = gemca;
    return 0;

fail:
    osh_gemca_prepared_free(gemca);
    return -1;
}

void osh_geometry_workspace_free(struct osh_geometry_workspace *ws) {
    size_t i;

    if (!ws) {
        return;
    }

    if (ws->bodies) {
        for (i = 0u; i < ws->nbodies; ++i) {
            geometry_body_free_fields(&ws->bodies[i]);
        }
        free(ws->bodies);
    }

    if (ws->zones) {
        for (i = 0u; i < ws->nzones; ++i) {
            geometry_zone_free_fields(&ws->zones[i]);
        }
        free(ws->zones);
    }

    if (ws->prepared) {
        osh_gemca_prepared_free(ws->prepared);
    }

    free(ws);
}
