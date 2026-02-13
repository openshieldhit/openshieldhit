#ifndef OSH_GEMCA2_CALC_SURFACE_H
#define OSH_GEMCA2_CALC_SURFACE_H

#include "gemca/osh_gemca2.h"

int osh_gemca2_add_surfaces(struct body *b, int n);
int osh_gemca2_add_surf_pars(struct surface *s, int type);
int osh_gemca2_check_surface_side(struct surface const *sf, struct ray const *r);

#endif /* OSH_GEMCA2_CALC_SURFACE_H */