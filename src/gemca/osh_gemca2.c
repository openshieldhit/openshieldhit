#include "gemca/osh_gemca2.h"

#include <stdio.h>
#include <stdlib.h>

#include "common/osh_logger.h"
#include "gemca/osh_gemca2_calc_body.h"
#include "gemca/osh_gemca2_calc_zone.h"
#include "gemca/osh_gemca2_defines.h"
#include "gemca/osh_gemca2_dist.h"
#include "gemca/parse/osh_gemca2_parse.h"

int osh_gemca_workspace_init(struct gemca_workspace **wg) {

    *wg = malloc(sizeof(struct gemca_workspace));
    if (*wg == NULL) {
        osh_alloc_failed("osh_gemca_workspace_init()");
        return 0;
    }
    return 1;
}

int osh_gemca_workspace_free(struct gemca_workspace *wg) {

    size_t i;

    for (i = 0; i < wg->nbodies; i++) {
        free(wg->bodies[i]->name);
        free(wg->bodies[i]->a);
        free(wg->bodies[i]);
    }
    free(wg->bodies);

    for (i = 0; i < wg->nzones; i++) {
        free(wg->zones[i]->name);
        free(wg->zones[i]);
    }
    free(wg->zones);

    free(wg);
    return 0;
}

int osh_gemca_load(const char *filename, struct gemca_workspace *g) {

    osh_gemca_parse(filename, g);

    printf("--- SETUP BODIES \n");
    osh_gemca_body_setup(g);
    printf("--- SETUP BODIES COMPLETED ---- \n\n");

    return 1;
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

    printf("Gemca->nBodies %lli ->nZones%lli\n", (unsigned long long) g->nbodies, (unsigned long long) g->nzones);

    for (i = 0; i < g->nbodies; i++) {
        osh_gemca_print_body(g->bodies[i]);
    }

    printf("\n");
    for (i = 0; i < g->nzones; i++) {
        osh_gemca_print_zone(g->zones[i]);
    }
    printf("\n");
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

    printf("----- PRINT BODY ----------------------------\n");
    printf("    Body name   : '%s'\n", b->name);
    printf("    Body type   : %i\n", b->type);
    printf("    Body nargs  : %i\n", b->na);
    printf("    Body nsurfs : %i\n", b->nsurfs);
    printf("    Body args   : ");
    for (i = 0; i < b->na; i++) {
        printf("%.2f ", b->a[i]);
    }
    printf("\n");
    printf("    Body surfaces... \n");
    for (i = 0; i < b->nsurfs; i++) {
        printf("         Surface %i is of type %i ", i, b->surfs[i]->type);
        printf("   parameters: ");
        for (j = 0; j < b->surfs[i]->np; j++) {
            printf("%.2f ", b->surfs[i]->p[j]);
        }
        printf("\n");
    }

    printf("\n");
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
    printf("----- PRINT ZONE ----------------------------\n");
    printf("    Zone name   : '%s'\n", z->name);
    printf("    Zone id     :  %llu\n", (unsigned long long) z->id);
    printf("    Zone medium :  %llu\n", (unsigned long long) z->medium);
    printf("    Zone tree follows...\n");
    osh_gemca_print_cgnodes(&z->node);
}

void osh_gemca_print_surface(struct surface const *s) {
    int i;

    printf("----- PRINT SURFACE -------------------------\n");
    printf("    Surface type : %i\n", s->type);
    printf("    Surface np   : %i\n", s->np);
    printf("    Surface params: ");
    for (i = 0; i < s->np; i++) {
        printf("%.2f ", s->p[i]);
    }
    printf("\n\n");
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
    printf("        This node pointer   : %p\n", (void *) self);

    if (self->type == _OSH_GEMCA_CGNODE_BODY) {
        printf("        Node type           : BODY '%s'\n", self->body->name);
        // printf("   '%s'\n", self->body->name);
    } else {
        printf("        Node type           : CGNODE\n");
        printf("        ->Left * -  Right    : %p '%c' %p\n", (void *) self->left, self->op, (void *) self->right);
    }
    printf("\n");

    if (self->left != NULL) {
        osh_gemca_print_cgnodes(self->left);
    }
    if (self->right != NULL) {
        osh_gemca_print_cgnodes(self->right);
    }
}
