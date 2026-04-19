#include "gemca/osh_gemca2.h"

#include <stdio.h>
#include <stdlib.h>

#include "common/osh_diag.h"
#include "gemca/osh_gemca2_calc_body.h"
#include "gemca/osh_gemca2_calc_zone.h"
#include "gemca/osh_gemca2_defines.h"
#include "gemca/osh_gemca2_dist.h"
#include "gemca/osh_gemca2_internal.h"

enum osh_status osh_gemca_prepared_free(struct osh_gemca_prepared *wg) {

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

/** @brief Print all bodies and zones in the workspace to the debug log. */
void osh_gemca_prepared_print(struct osh_gemca_prepared const *g, struct osh_diag_sink const *diag) {

    size_t i;

    if (!g || !diag || !diag->emit || diag->min_level > OSH_DIAG_LEVEL_DEBUG) {
        return;
    }

    OSH_DIAG_DEBUGF(
        diag, "Gemca->nBodies %llu ->nZones %llu", (unsigned long long) g->nbodies, (unsigned long long) g->nzones);

    for (i = 0; i < g->nbodies; i++) {
        osh_gemca_print_body(g->bodies[i], diag);
    }

    for (i = 0; i < g->nzones; i++) {
        osh_gemca_print_zone(g->zones[i], diag);
    }
}

/** @brief Print body parameters (name, type, args, surfaces) to the debug log. */
void osh_gemca_print_body(struct body const *b, struct osh_diag_sink const *diag) {

    int i;
    int j;

    if (!b || !diag || !diag->emit || diag->min_level > OSH_DIAG_LEVEL_DEBUG) {
        return;
    }

    OSH_DIAG_DEBUGF(diag, "%s", "------------------------------------------------------------");
    OSH_DIAG_DEBUGF(diag, "    Body name   : '%s'", b->name);
    OSH_DIAG_DEBUGF(diag, "    Body type   : %i", b->type);
    OSH_DIAG_DEBUGF(diag, "    Body nargs  : %i", b->na);
    OSH_DIAG_DEBUGF(diag, "    Body nsurfs : %i", b->nsurfs);
    OSH_DIAG_DEBUGF(diag, "    Body args   : ");
    for (i = 0; i < b->na; i++) {
        OSH_DIAG_DEBUGF(diag, "%.2f ", b->a[i]);
    }
    OSH_DIAG_DEBUGF(diag, "    Body surfaces...");
    for (i = 0; i < b->nsurfs; i++) {
        OSH_DIAG_DEBUGF(diag, "        Surface %i type %i  params:", i, b->surfs[i]->type);
        for (j = 0; j < b->surfs[i]->np; j++) {
            OSH_DIAG_DEBUGF(diag, " %.2f", b->surfs[i]->p[j]);
        }
    }
}

/** @brief Print zone name, material, and CSG tree to the debug log. */
void osh_gemca_print_zone(struct zone const *z, struct osh_diag_sink const *diag) {
    if (!z || !diag || !diag->emit || diag->min_level > OSH_DIAG_LEVEL_DEBUG) {
        return;
    }

    OSH_DIAG_DEBUGF(diag, "%s", "------------------------------------------------------------");
    OSH_DIAG_DEBUGF(diag, "    Zone name   : '%s'", z->name);
    OSH_DIAG_DEBUGF(diag, "    Zone material :  %s", z->material_name ? z->material_name : "(unset)");
    OSH_DIAG_DEBUGF(diag, "    Zone material index :  %llu", (unsigned long long) z->material_idx);
    OSH_DIAG_DEBUGF(diag, "    Zone tree follows...");
    osh_gemca_print_cgnodes(&z->node, diag);
}

void osh_gemca_print_surface(struct surface const *s, struct osh_diag_sink const *diag) {
    int i;

    if (!s || !diag || !diag->emit || diag->min_level > OSH_DIAG_LEVEL_DEBUG) {
        return;
    }

    OSH_DIAG_DEBUGF(diag, "%s", "------------------------------------------------------------");
    OSH_DIAG_DEBUGF(diag, "    Surface type : %i", s->type);
    OSH_DIAG_DEBUGF(diag, "    Surface np   : %i", s->np);
    OSH_DIAG_DEBUGF(diag, "    Surface params:");
    for (i = 0; i < s->np; i++) {
        OSH_DIAG_DEBUGF(diag, " %.2f", s->p[i]);
    }
    OSH_DIAG_DEBUGF(diag, "%s", "");
}

/** @brief Recursively print the CSG node tree rooted at `self` to the debug log. */
void osh_gemca_print_cgnodes(struct cgnode const *self, struct osh_diag_sink const *diag) {
    if (!self || !diag || !diag->emit || diag->min_level > OSH_DIAG_LEVEL_DEBUG) {
        return;
    }

    OSH_DIAG_DEBUGF(diag, "        This node pointer   : %p", (void *) self);

    if (self->type == _OSH_GEMCA_CGNODE_BODY) {
        OSH_DIAG_DEBUGF(diag, "        Node type           : BODY '%s'", self->body->name);
    } else {
        OSH_DIAG_DEBUGF(diag, "        Node type           : CGNODE");
        OSH_DIAG_DEBUGF(
            diag, "        ->Left * -  Right    : %p '%c' %p", (void *) self->left, self->op, (void *) self->right);
    }
    OSH_DIAG_DEBUGF(diag, "%s", "");

    if (self->left != NULL) {
        osh_gemca_print_cgnodes(self->left, diag);
    }
    if (self->right != NULL) {
        osh_gemca_print_cgnodes(self->right, diag);
    }
}
