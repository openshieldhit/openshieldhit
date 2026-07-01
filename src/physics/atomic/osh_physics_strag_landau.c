#include "physics/atomic/osh_physics_strag_landau.h"

#include <math.h>
#include <stddef.h>

#include "physics/atomic/osh_physics_strag_landau_coeffs.h"

/* Per-band transform x'(u); ids match emit_landau.py (0=-log u, 1=u, 2=-log(1-u)). */
static double _lan_xform(double u, int kind) {
    switch (kind) {
    case 0:
        return -log(u);
    case 1:
        return u;
    default:
        return -log1p(-u);
    }
}

double osh_physics_strag_landau_lambda(double u) {
    int ub;          /* selected u-band index (0..OSH_LAN_NB-1)           */
    int k;           /* polynomial power in the Horner loop               */
    double xn;       /* transformed u mapped to ~[-1,1] within the u-band */
    double lam;      /* result: reduced (universal) Landau variable λ     */
    double const *c; /* coefficient row for the selected u-band           */

    if (u < OSH_LAN_UMIN) {
        u = OSH_LAN_UMIN;
    }
    if (u > OSH_LAN_UMAX) {
        u = OSH_LAN_UMAX;
    }

    for (ub = 0; ub < OSH_LAN_NB - 1; ++ub) {
        if (u <= osh_lan_uhi[ub]) {
            break;
        }
    }
    xn = (_lan_xform(u, osh_lan_kind[ub]) - osh_lan_xmid[ub]) / osh_lan_xhalf[ub];

    /* λ = Horner in xn over powers 0..deg. */
    c = &osh_lan_coef[(size_t) ub * OSH_LAN_MAXP];
    lam = 0.0;
    for (k = osh_lan_deg[ub]; k >= 0; --k) {
        lam = lam * xn + c[k];
    }
    return lam;
}
