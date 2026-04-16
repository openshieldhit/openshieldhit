#ifndef OSH_APP_OSH_GEOMETRY_PARSE_KEYS_H
#define OSH_APP_OSH_GEOMETRY_PARSE_KEYS_H

/**
 * @file osh_geometry_parse_keys.h
 * @brief OpenShieldHIT `geo.dat` file-syntax keys used by the app parser.
 *
 * @details
 * These are app-private parsing tokens. They intentionally do not live in
 * `osh_core`, because file syntax belongs to `src/apps/osh/`, not to the
 * public geometry API or the GEMCA prepare/runtime layers.
 */

#define OSH_GEO_KEY_END "end"
#define OSH_GEO_KEY_ASSIGNMA "assignma"
#define OSH_GEO_KEY_ASSIGNMAT "assignmat"

#endif /* OSH_APP_OSH_GEOMETRY_PARSE_KEYS_H */
