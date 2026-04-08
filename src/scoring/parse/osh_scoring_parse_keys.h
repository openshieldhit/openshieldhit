#ifndef OSH_SCORING_PARSE_KEYS_H
#define OSH_SCORING_PARSE_KEYS_H

/* clang-format off */

/* ---- Filter section keys ------------------------------------------------- */
#define OSH_SCORING_KEY_NAME        "name"

/* Filter rule fields (case-insensitive, stored lower-case in dispatch table) */
#define OSH_SCORING_KEY_FILTER_Z    "z"
#define OSH_SCORING_KEY_FILTER_A    "a"
#define OSH_SCORING_KEY_FILTER_E    "e"
#define OSH_SCORING_KEY_FILTER_GEN  "gen"
#define OSH_SCORING_KEY_FILTER_ID   "id"

/* ---- Settings section keys ----------------------------------------------- */
/* OSH_SCORING_KEY_NAME shared */
#define OSH_SCORING_KEY_RESCALE     "rescale"
#define OSH_SCORING_KEY_OFFSET      "offset"
#define OSH_SCORING_KEY_MEDIUM      "medium"
#define OSH_SCORING_KEY_NKMEDIUM    "nkmedium"
#define OSH_SCORING_KEY_SITEDIAM    "sitediameter"  /* canonical; "sitediam" is alias */
#define OSH_SCORING_KEY_DENSITY     "density"       /* canonical; "rho" is alias */
#define OSH_SCORING_KEY_MAXCOUNT    "maxcount"      /* canonical; "npart" is alias */

/* ---- Geometry section keys ----------------------------------------------- */
/* OSH_SCORING_KEY_NAME shared */
#define OSH_SCORING_KEY_GEO_X       "x"
#define OSH_SCORING_KEY_GEO_Y       "y"
#define OSH_SCORING_KEY_GEO_Z       "z"
#define OSH_SCORING_KEY_GEO_R       "r"
#define OSH_SCORING_KEY_GEO_ROT     "rotation"      /* canonical; "rot" is alias */
#define OSH_SCORING_KEY_GEO_ZONES   "zones"

/* ---- Output section keys ------------------------------------------------- */
#define OSH_SCORING_KEY_FILENAME    "filename"
#define OSH_SCORING_KEY_GEO_REF     "geo"
#define OSH_SCORING_KEY_FILEFORMAT  "fileformat"    /* canonical; "format" is alias */
#define OSH_SCORING_KEY_QUANTITY    "quantity"

/* clang-format on */

#endif /* OSH_SCORING_PARSE_KEYS_H */
