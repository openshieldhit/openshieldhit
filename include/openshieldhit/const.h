#ifndef OPENSHIELDHIT_CONST_H
#define OPENSHIELDHIT_CONST_H

/**
 * @file openshieldhit/const.h
 * @brief Physical and mathematical constants used throughout OpenShieldHIT.
 *
 * @details
 * All values are taken from NIST (https://physics.nist.gov/cuu/Constants/).
 * Notable 2020 redefinitions:
 *   OSH_AMU       was 931.4940954  MeV/c²
 *   OSH_NAVOGADRO was 6.022140857e23 mol⁻¹
 */

#define OSH_AMU 931.49410242         /* 1 Da (dalton, unified atomic mass unit) [MeV/c²] */
#define OSH_NAVOGADRO 6.02214076e23  /* Avogadro constant [mol⁻¹] */
#define OSH_MEVG2GY 1.602176634e-10  /* 1 MeV/g in J/kg */
#define OSH_EV2JOULE 1.602176634e-19 /* 1 eV in J */

#define OSH_M_PI 3.14159265358979323846       /* pi */
#define OSH_M_1_PI 0.31830988618379067154     /* 1/pi */
#define OSH_M_PI_180 0.0174532925199432957692 /* pi/180 */

#define OSH_SIGMA2FWHM 2.3548200450309493  /* 1σ → FWHM of Gaussian: 2·sqrt(2·ln(2)) */
#define OSH_FWHM2SIGMA 0.42466090014400953 /* 1/OSH_SIGMA2FWHM */

#endif /* OPENSHIELDHIT_CONST_H */
