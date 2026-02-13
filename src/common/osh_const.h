/* Various physical and mathematical constants */
#ifndef _OSH_CONST
#define _OSH_CONST

/* taken from NIST: https://physics.nist.gov/cuu/Constants/index.html */
/* OSH_AMU         Changed in 2020, old definition was 931.4940954  MeV/c**2 */
/* OSH_NAVOGADRO   Changed in 2020, old definition was 6.022140857e23  mol^-1 */

#define OSH_AMU 931.49410242         /* 1 amu in [MeV/c**2] */
#define OSH_NAVOGADRO 6.02214076e23  /* Avogadro constant in [mol^-1] */
#define OSH_MEVG2GY 1.602176634e-10  /* 1 MeV/g in J/kg */
#define OSH_EV2JOULE 1.602176634e-19 /* 1 eV in J */

#define OSH_M_PI 3.14159265358979323846       /* pi (WolframAlpha) */
#define OSH_M_1_PI 0.31830988618379067154     /* 1/pi (WolframAlpha) */
#define OSH_M_PI_180 0.0174532925199432957692 /* pi/180 (WolframAlpha) */

#define OSH_SIGMA2FWHM 2.3548200450309493  /* factor to convert 1 sigma to FWHM of Gaussian dist 2*sqrt(2*ln(2)) */
#define OSH_FWHM2SIGMA 0.42466090014400953 /* = 1.0 / OSH_SIGMA2FWHM */
#endif                                     /* !_OSH_CONST */
