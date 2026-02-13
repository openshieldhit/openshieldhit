#include "gemca/voxel/osh_gemca2_voxel_hu.h"

#include "common/osh_interpolate.h"
#include "gemca/voxel/osh_gemca2_voxel_defines.h"
#include "gemca/voxel/osh_voxel_mat_schneider2000.h"

static inline float _wepl_minohara1993(int16_t hu);
static inline float _wepl_jacob1996(int16_t hu);
static inline float _wepl_geiss1999(int16_t hu);

float osh_gemca_voxel_hu2rho(int16_t hu, char alg) {
    float ret = 0.0;

    // TODO:  leszek did some strange test if HU >= 10000 ? Looks like a mode switch
    (void) (alg); // unused

    if (hu < -1000) {
        hu = -1000;
    }
    if (hu > 1600) {
        hu = 1600;
    }

    if (hu <= 98) {
        ret = 1.03091 + 0.0010297 * hu;
    } else if (hu <= 14) {
        ret = 1.018 + 0.893 * 0.001 * hu;
    } else if (hu <= 23) {
        ret = 1.03;
    } else if (hu <= 100) {
        ret = 1.003 + 1.169 * 0.001 * hu;
    } else {
        ret = 1.017 + 0.592 * 0.001 * hu;
    }
    return ret;
}

/* lookup material index from HU value */
/* hu must be between -1000 and 1600 */
int osh_gemca_voxel_hu2idx(int16_t hu) {
    if (hu < -1000) {
        return 0.0;
    }
    if (hu > 1600) {
        return 0.0;
    }
    return osh_binary_search_i2(hu, _ct_hu, _nmat + 1);
}

float osh_gemca_voxel_hu2wepl(int16_t hu, char alg) {

    float wepl = 0.0;

    if (hu < -1000) {
        return 0.0;
    }

    if (hu > 4000) {
        return 0.0;
    }

    switch (alg) {
    case OSH_GEMCA_VOXEL_HU2WEPL_ALG1:
        wepl = _wepl_minohara1993(hu);
        break;
    case OSH_GEMCA_VOXEL_HU2WEPL_ALG2:
        wepl = _wepl_jacob1996(hu);
        break;
    case OSH_GEMCA_VOXEL_HU2WEPL_ALG3:
        wepl = _wepl_geiss1999(hu);
        break;
    default:
        return 0.0;
        break;
    }
    return wepl * 1000.0;
}

static inline float _wepl_minohara1993(int16_t hu) {
    if (hu < -49) {
        return 1.075e-3 * hu + 1.050;
    } else {
        return 4.597e-4 * hu + 1.019;
    }
}

static inline float _wepl_jacob1996(int16_t hu) {
    if (hu < -60.81) { /* float to int comparison */
        return 1.011e-3 * hu + 1.052;
    } else {
        return 4.190e-4 * hu + 1.016;
    }
}

static inline float _wepl_geiss1999(int16_t hu) {
    // TODO
    (void) (hu); // unused
    return 0.0;
}
