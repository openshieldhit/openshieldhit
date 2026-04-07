#include "gemca/osh_gemca2.h"

#include <stdio.h>
#include <stdlib.h>

#include "common/osh_logger.h"
#include "gemca/osh_gemca2_calc_body.h"
#include "gemca/osh_gemca2_calc_zone.h"
#include "gemca/osh_gemca2_defines.h"
#include "gemca/osh_gemca2_dist.h"
#include "gemca/parse/osh_gemca2_parse.h"

enum osh_status osh_gemca_workspace_init(struct gemca_workspace **wg) {

    *wg = calloc(1, sizeof(struct gemca_workspace));
    if (*wg == NULL) {
        osh_alloc_failed("osh_gemca_workspace_init()");
        return OSH_ENOMEM;
    }
    return OSH_OK;
}

enum osh_status osh_gemca_workspace_free(struct gemca_workspace *wg) {

    size_t i;

    if (wg == NULL) {
        return OSH_OK;
    }

    if (wg->bodies != NULL) {
        for (i = 0; i < wg->nbodies; i++) {
            if (wg->bodies[i] != NULL) {
                free(wg->bodies[i]->name);
                free(wg->bodies[i]->a);
                free(wg->bodies[i]);
            }
        }
    }
    free((void *) wg->bodies);

    if (wg->zones != NULL) {
        for (i = 0; i < wg->nzones; i++) {
            if (wg->zones[i] != NULL) {
                free(wg->zones[i]->name);
                free(wg->zones[i]->material_name);
                free(wg->zones[i]);
            }
        }
    }
    free((void *) wg->zones);

    free(wg->filename);
    free(wg);
    return OSH_OK;
}

enum osh_status osh_gemca_load(char const *filename, struct gemca_workspace *g) {
    enum osh_status rc;

    rc = osh_gemca_parse(filename, g);
    if (rc != OSH_OK) {
        return rc;
    }

    osh_debug("Setting up bodies...");
    rc = osh_gemca_body_setup(g);
    if (rc != OSH_OK) {
        return rc;
    }
    osh_debug("Body setup complete.");

    return OSH_OK;
}

/** @brief Print all bodies and zones in the workspace to the debug log. */
void osh_gemca_print_gemca(struct gemca_workspace const *g) {

    size_t i;

    osh_debug("Gemca->nBodies %llu ->nZones %llu", (unsigned long long) g->nbodies, (unsigned long long) g->nzones);

    for (i = 0; i < g->nbodies; i++) {
        osh_gemca_print_body(g->bodies[i]);
    }

    for (i = 0; i < g->nzones; i++) {
        osh_gemca_print_zone(g->zones[i]);
    }
}

/** @brief Print body parameters (name, type, args, surfaces) to the debug log. */
void osh_gemca_print_body(struct body const *b) {

    int i;
    int j;

    osh_debug(OSH_LOG_HLINE);
    osh_debug("    Body name   : '%s'", b->name);
    osh_debug("    Body type   : %i", b->type);
    osh_debug("    Body nargs  : %i", b->na);
    osh_debug("    Body nsurfs : %i", b->nsurfs);
    osh_debug("    Body args   : ");
    for (i = 0; i < b->na; i++) {
        osh_debug("%.2f ", b->a[i]);
    }
    osh_debug("    Body surfaces...");
    for (i = 0; i < b->nsurfs; i++) {
        osh_debug("    " OSH_LOG_INDENT "Surface %i type %i  params:", i, b->surfs[i]->type);
        for (j = 0; j < b->surfs[i]->np; j++) {
            osh_debug(" %.2f", b->surfs[i]->p[j]);
        }
    }
}

/** @brief Print zone name, material, and CSG tree to the debug log. */
void osh_gemca_print_zone(struct zone const *z) {
    osh_debug(OSH_LOG_HLINE);
    osh_debug("    Zone name   : '%s'", z->name);
    osh_debug("    Zone material :  %s", z->material_name ? z->material_name : "(unset)");
    osh_debug("    Zone material index :  %llu", (unsigned long long) z->material_idx);
    osh_debug("    Zone tree follows...");
    osh_gemca_print_cgnodes(&z->node);
}

void osh_gemca_print_surface(struct surface const *s) {
    int i;

    osh_debug(OSH_LOG_HLINE);
    osh_debug("    Surface type : %i", s->type);
    osh_debug("    Surface np   : %i", s->np);
    osh_debug("    Surface params:");
    for (i = 0; i < s->np; i++) {
        osh_debug(" %.2f", s->p[i]);
    }
    osh_debug("%s", "");
}

/** @brief Recursively print the CSG node tree rooted at `self` to the debug log. */
void osh_gemca_print_cgnodes(struct cgnode const *self) {
    osh_debug("        This node pointer   : %p", (void *) self);

    if (self->type == _OSH_GEMCA_CGNODE_BODY) {
        osh_debug("        Node type           : BODY '%s'", self->body->name);
    } else {
        osh_debug("        Node type           : CGNODE");
        osh_debug("        ->Left * -  Right    : %p '%c' %p", (void *) self->left, self->op, (void *) self->right);
    }
    osh_debug("%s", "");

    if (self->left != NULL) {
        osh_gemca_print_cgnodes(self->left);
    }
    if (self->right != NULL) {
        osh_gemca_print_cgnodes(self->right);
    }
}
