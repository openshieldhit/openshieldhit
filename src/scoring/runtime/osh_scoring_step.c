#include "scoring/runtime/osh_scoring_step.h"

#include <stdlib.h>
#include <string.h>

#include "common/raytrace/osh_raytrace.h"

static int axis_index(struct osh_scoring_geometry_runtime const *geo, char const *label);
static enum osh_status mesh_geometry_to_grid(struct osh_scoring_geometry_runtime const *geo,
                                             struct osh_raytrace_grid *grid,
                                             double *voxel_volume_out);
static enum osh_status score_group_energy(struct osh_scoring_runtime *rt,
                                          struct osh_scoring_geometry_score_group const *group,
                                          struct osh_voxel_crossing const *crossings,
                                          size_t ncross,
                                          struct particle const *part,
                                          struct step const *st);
static enum osh_status score_group_fluence(struct osh_scoring_runtime *rt,
                                           struct osh_scoring_geometry_score_group const *group,
                                           struct osh_voxel_crossing const *crossings,
                                           size_t ncross,
                                           struct particle const *part,
                                           struct step const *st,
                                           double voxel_volume);

enum osh_status
osh_scoring_score_step(struct osh_scoring_runtime *rt, struct particle const *part, struct step const *st) {
    size_t i;
    size_t cap;
    size_t ncross;
    double voxel_volume;
    struct osh_raytrace_grid grid;
    struct osh_voxel_crossing *crossings;
    enum osh_status rc;
    int hit;

    if (!rt || !part || !st) {
        return OSH_EINVAL;
    }
    if (!(st->ds > 0.0)) {
        return OSH_EINVAL;
    }

    for (i = 0; i < rt->ngeometries; ++i) {
        struct osh_scoring_geometry_runtime const *geo = &rt->geometries[i];
        size_t g;

        rc = mesh_geometry_to_grid(geo, &grid, &voxel_volume);
        if (rc == OSH_ENOTSUP && geo->npages == 0u) {
            continue;
        }
        if (rc != OSH_OK) {
            return rc;
        }

        cap = grid.n[0] + grid.n[1] + grid.n[2];
        if (cap == 0u) {
            continue;
        }
        crossings = (struct osh_voxel_crossing *) calloc(cap, sizeof(*crossings));
        if (!crossings) {
            return OSH_ENOMEM;
        }

        hit = osh_raytrace_traverse(&grid, st->p, st->v, st->ds, crossings, &ncross);
        if (!hit || ncross == 0u) {
            free(crossings);
            continue;
        }

        /* TODO: hoist more geometry/crossing-derived work out of the per-group
         * loops: store inv_voxel_volume in geometry runtime, precompute
         * crossing fractions once per step, and later reuse scratch buffers
         * instead of allocating crossings[] on every call. */
        for (g = 0; g < geo->ngroups; ++g) {
            switch (geo->groups[g].score_kind) {
            case OSH_SCORING_SCORE_ENERGY:
                rc = score_group_energy(rt, &geo->groups[g], crossings, ncross, part, st);
                break;
            case OSH_SCORING_SCORE_FLUENCE:
                rc = score_group_fluence(rt, &geo->groups[g], crossings, ncross, part, st, voxel_volume);
                break;
            default:
                rc = OSH_ENOTSUP;
                break;
            }
            if (rc != OSH_OK) {
                free(crossings);
                return rc;
            }
        }

        free(crossings);
    }

    return OSH_OK;
}

enum osh_status
osh_scoring_score_point(struct osh_scoring_runtime *rt, struct particle const *part, struct position const *pos) {
    (void) rt;
    (void) part;
    (void) pos;
    return OSH_ENOTSUP;
}

static int axis_index(struct osh_scoring_geometry_runtime const *geo, char const *label) {
    size_t i;

    for (i = 0; i < geo->naxes; ++i) {
        if (strcmp(geo->axes[i].label, label) == 0) {
            return (int) i;
        }
    }
    return -1;
}

static enum osh_status mesh_geometry_to_grid(struct osh_scoring_geometry_runtime const *geo,
                                             struct osh_raytrace_grid *grid,
                                             double *voxel_volume_out) {
    int ix;
    int iy;
    int iz;
    double dx;
    double dy;
    double dz;

    if (!geo || !grid || !voxel_volume_out) {
        return OSH_EINVAL;
    }
    if (geo->geo_kind != OSH_SCORING_GEO_MESH) {
        return OSH_ENOTSUP;
    }
    if (geo->has_rotation) {
        return OSH_ENOTSUP;
    }
    if (geo->naxes != 3u) {
        return OSH_EINVAL;
    }

    ix = axis_index(geo, "X");
    iy = axis_index(geo, "Y");
    iz = axis_index(geo, "Z");
    if (ix < 0 || iy < 0 || iz < 0) {
        return OSH_EINVAL;
    }
    if (geo->axes[ix].nbins <= 0 || geo->axes[iy].nbins <= 0 || geo->axes[iz].nbins <= 0) {
        return OSH_EINVAL;
    }

    dx = (geo->axes[ix].hi - geo->axes[ix].lo) / (double) geo->axes[ix].nbins;
    dy = (geo->axes[iy].hi - geo->axes[iy].lo) / (double) geo->axes[iy].nbins;
    dz = (geo->axes[iz].hi - geo->axes[iz].lo) / (double) geo->axes[iz].nbins;
    if (!(dx > 0.0) || !(dy > 0.0) || !(dz > 0.0)) {
        return OSH_EINVAL;
    }

    grid->origin[0] = geo->axes[ix].lo;
    grid->origin[1] = geo->axes[iy].lo;
    grid->origin[2] = geo->axes[iz].lo;
    grid->spacing[0] = dx;
    grid->spacing[1] = dy;
    grid->spacing[2] = dz;
    grid->n[0] = (size_t) geo->axes[ix].nbins;
    grid->n[1] = (size_t) geo->axes[iy].nbins;
    grid->n[2] = (size_t) geo->axes[iz].nbins;
    *voxel_volume_out = dx * dy * dz;
    return OSH_OK;
}

static enum osh_status score_group_energy(struct osh_scoring_runtime *rt,
                                          struct osh_scoring_geometry_score_group const *group,
                                          struct osh_voxel_crossing const *crossings,
                                          size_t ncross,
                                          struct particle const *part,
                                          struct step const *st) {
    size_t i;
    size_t j;
    double frac;
    struct osh_scoring_page_runtime *page;

    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        if (!osh_scoring_page_passes_filters(rt, page, part, st)) {
            continue;
        }
        for (j = 0; j < ncross; ++j) {
            if (crossings[j].idx >= page->len) {
                return OSH_ESTATE;
            }
            frac = crossings[j].path_len / st->ds;
            page->data[crossings[j].idx] += st->de * frac;
        }
    }
    return OSH_OK;
}

static enum osh_status score_group_fluence(struct osh_scoring_runtime *rt,
                                           struct osh_scoring_geometry_score_group const *group,
                                           struct osh_voxel_crossing const *crossings,
                                           size_t ncross,
                                           struct particle const *part,
                                           struct step const *st,
                                           double voxel_volume) {
    size_t i;
    size_t j;
    struct osh_scoring_page_runtime *page;

    (void) st;

    for (i = 0; i < group->npages; ++i) {
        page = &rt->pages[group->first_page + i];
        if (!osh_scoring_page_passes_filters(rt, page, part, st)) {
            continue;
        }
        for (j = 0; j < ncross; ++j) {
            if (crossings[j].idx >= page->len) {
                return OSH_ESTATE;
            }
            page->data[crossings[j].idx] += crossings[j].path_len / voxel_volume;
        }
    }
    return OSH_OK;
}
