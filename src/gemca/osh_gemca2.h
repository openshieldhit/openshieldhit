#ifndef _OSH_GEMCA2
#define _OSH_GEMCA2

#include <float.h>
#include <stdio.h>

#include "transport/osh_transport.h" /* needed for struct ray */

#define _OSH_GEMCA_CGNODE_BODY 0
#define _OSH_GEMCA_CGNODE_COMPOSITE 1

#ifdef INFINITY
#define OSH_GEMCA_INFINITY INFINITY
#else
#define OSH_GEMCA_INFINITY 1e300
#endif

/* can be changed if needed */
#define OSH_GEMCA_SMALL 1e-12
#define OSH_GEMCA_STEPLIM                                                                                              \
    1e-8 /* minimal step to avoid getting stuck on                                                                     \
               surface due to numerical precision */

struct gemca_workspace {  /* workspace for gemca */
    struct body **bodies; /* list of pointers to all bodies */
    struct zone **zones;  /* list of pointers to zones */
    size_t nbodies;       /* total number of bodies */
    size_t nzones;        /* total number of zones */
    char *filename;       /* path to the geo.dat file */
};

struct body {               /* a body primitive */
    double t[16];           /* 4x4 transformation matrix for translating OSH_COORD_UNIVERSE
                               --> OSH_COORD_B**** */
    struct surface **surfs; /* pointer to list of surfaces */
    size_t lineno;          /* line number where this body was defined */
    char *name;             /* user given name for this body */
    char *filename_vox;     /* placeholder for path to voxel file, in case of vox body
                             */
    double *a;              /* list of arguments given to this body */
    int na;                 /* number of arguments in list *a */
    int type;               /* type identifier */
    int nsurfs;             /* number of surfaces */
    char coord;             /* body parameters are in this coordinate system */
};

struct cgnode {
    double bb_max[3];     /* bounding box max // TODO  */
    double bb_min[3];     /* bounding box min // TODO */
    struct cgnode *left;  /* used only if this is a composite node */
    struct cgnode *right; /* used only if this is a composite node */
    struct body *body;    /* used only if this is a leaf (=body) node */
    int type;             /* _OSH_GEMCA_CGNODE_* : mark if this is a leaf (=body) node or a
                             composite node */
    char op;              /* operator, if this is a composite node */
    char _is_inside;      /* updated while raycasting to see whether ray is inside or
                             outside node */
};

struct zone {           /* zone description */
    struct cgnode node; /* top level node of the abstrac syntax tree which holds
                           the body description */
    size_t id;          /* number of this zone, starting at 1 */
    size_t lineno;      /* fist line number where this zone was defined */
    size_t medium;      /* medium/material ID of this zone */
    size_t ntokens;     /* number of tokens */
    char **tokens;      /* list of tokens */
    char *name;         /* user given name of this zone */
};

struct surface {  /* surface descriptions */
    double *p;    /* list of parameters describing the surface, what goes in here
                     depends on the type */
    double _dist; /* calculated distance to zone, negative if ray does not cross */
    int np;       /* number of parameters in list *p */
    int type;     /* type identifier */
};

/* loads filename and prepares *gemca workspace */
int osh_gemca_workspace_init(struct gemca_workspace **wg);
int osh_gemca_workspace_free(struct gemca_workspace *wg);
int osh_gemca_load(const char *filename, struct gemca_workspace *g);

/* for a given ray and *g workspace, return what zone ID (starts at 1) we are in
 */
size_t osh_gemca_zone(struct gemca_workspace g, struct ray r);
size_t osh_gemca_zone_index(struct gemca_workspace g, struct ray r); /* return index of zone, can be used
                                                                       directly on g->zones[index]*/

/* for a given ray and *g workspace, return distance to nearest surface along
 * ray */
double osh_gemca_dist(struct zone *z, struct ray const *r);

void osh_gemca_print_gemca(struct gemca_workspace const *g); /* print entire workspace */
void osh_gemca_print_body(struct body const *b);
void osh_gemca_print_zone(struct zone const *z);
void osh_gemca_print_surface(struct surface const *s);
void osh_gemca_print_cgnodes(struct cgnode const *self);

#endif /* gemca2.h */
