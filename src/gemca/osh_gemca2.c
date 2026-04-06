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

    osh_debug("setting up bodies");
    rc = osh_gemca_body_setup(g);
    if (rc != OSH_OK) {
        return rc;
    }
    osh_debug("body setup complete");

    return OSH_OK;
}

/* for a given ray and *g workspace, return what zone we are in */
size_t osh_gemca_zone(struct gemca_workspace g, struct ray r) {
    size_t zone;

    zone = osh_gemca_get_zone(&g, &r);
    return zone;
}

/* for a given ray and *g workspace, return what zone we are in */
size_t osh_gemca_zone_index(struct gemca_workspace g, struct ray r) {
    size_t zidx;

    zidx = osh_gemca_get_zone_index(&g, &r);
    return zidx;
}

/* for a given ray and zone workspace, return distance to nearest surface along
 * ray */
double osh_gemca_dist(struct zone *z, struct ray const *r) {
    double d;

    d = osh_gemca_get_distance(z, r); // TODO: this double calling is just during debugging
    return d;
}

/**
 * @brief Print the gemca workspace
 *
 * @details
 *
 * @param[in]
 * @param[in]
 *
 * @returns
 *
 * @author Niels Bassler
 */
void osh_gemca_print_gemca(struct gemca_workspace const *g) {

    size_t i;

    osh_info("Gemca->nBodies %llu ->nZones%llu", (unsigned long long) g->nbodies, (unsigned long long) g->nzones);

    for (i = 0; i < g->nbodies; i++) {
        osh_gemca_print_body(g->bodies[i]);
    }

    osh_info("%s", "");
    for (i = 0; i < g->nzones; i++) {
        osh_gemca_print_zone(g->zones[i]);
    }
    osh_info("%s", "");
}

/**
 * @brief Print a given body
 *
 * @details
 *
 * @param[in]
 * @param[in]
 *
 * @returns
 *
 * @author Niels Bassler
 */
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
    osh_debug("%s", "");
    osh_debug("    Body surfaces...");
    for (i = 0; i < b->nsurfs; i++) {
        osh_debug("    " OSH_LOG_INDENT "Surface %i type %i  params:", i, b->surfs[i]->type);
        for (j = 0; j < b->surfs[i]->np; j++) {
            osh_debug(" %.2f", b->surfs[i]->p[j]);
        }
        osh_debug("%s", "");
    }
    osh_debug("%s", "");
}

/**
 * @brief Print the zone data.
 *
 * @details
 *
 * @param[in]
 * @param[in]
 *
 * @returns
 *
 * @author Niels Bassler
 */
void osh_gemca_print_zone(struct zone const *z) {
    osh_debug(OSH_LOG_HLINE);
    osh_debug("    Zone name   : '%s'", z->name);
    osh_debug("    Zone id     :  %llu", (unsigned long long) z->id);
    osh_debug("    Zone medium :  %llu", (unsigned long long) z->medium);
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

/**
 * @brief Print the abstract syntax tree.
 *
 * @details
 *
 * @param[in]
 * @param[in]
 *
 * @returns
 *
 * @author Niels Bassler
 */
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
