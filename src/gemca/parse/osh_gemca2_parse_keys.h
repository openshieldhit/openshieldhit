#ifndef _OSH_GEMCA2_KEYS
#define _OSH_GEMCA2_KEYS

/* list of keys used in geo.dat */
#define OSH_GEMCA_KEY_END "end"
#define OSH_GEMCA_KEY_SPH "sph"
#define OSH_GEMCA_KEY_WED "wed"
#define OSH_GEMCA_KEY_ARB "arb"
#define OSH_GEMCA_KEY_BOX "box"
#define OSH_GEMCA_KEY_VOX "vox"
#define OSH_GEMCA_KEY_RPP "rpp"
#define OSH_GEMCA_KEY_RCC "rcc"
#define OSH_GEMCA_KEY_REC "rec"
#define OSH_GEMCA_KEY_TRC "trc"
#define OSH_GEMCA_KEY_ELL "ell"

#define OSH_GEMCA_KEY_YZP "yzp"
#define OSH_GEMCA_KEY_XZP "xzp"
#define OSH_GEMCA_KEY_XYP "xyp"
#define OSH_GEMCA_KEY_PLA "pla"

#define OSH_GEMCA_KEY_ROT "rot"
#define OSH_GEMCA_KEY_CPY "cpy"
#define OSH_GEMCA_KEY_MOV "mov"

#define OSH_GEMCA_KEY_ASSIGNMA "assignma"   /* FLUKA compatible key for material assignment (no 't' at the end)*/
#define OSH_GEMCA_KEY_ASSIGNMAT "assignmat" /* allow both spellings for compatibility */

/* zone definition may have the "OR " in the arguments */
/* Currently these are not needed, since they are hardcoded for better readability. */
/*
 #define OSH_GEMCA_KEY_SL_OR      "or"  // notice this is a string, not a char, unlike the rest
 #define OSH_GEMCA_KEY_CL_INCL    '+'
 #define OSH_GEMCA_KEY_CL_EXCL    '-'
 #define OSH_GEMCA_KEY_CL_PIPE    '|'
 #define OSH_GEMCA_KEY_CL_PARO    '('
 #define OSH_GEMCA_KEY_CL_PARC    ')'
 */

#endif /* _OSH_GEMCA2_KEYS */
