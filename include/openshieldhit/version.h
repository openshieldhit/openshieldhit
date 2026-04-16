#ifndef OPENSHIELDHIT_VERSION_H
#define OPENSHIELDHIT_VERSION_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file version.h
 * @brief Public version query API for the osh_core library.
 */

/** @returns The full version string (e.g. "0.3.1-5-gabcdef-dirty"). */
char const *osh_version_string(void);

/** @returns The major version component. */
int osh_version_major(void);

/** @returns The minor version component. */
int osh_version_minor(void);

/** @returns The patch version component. */
int osh_version_patch(void);

#ifdef __cplusplus
}
#endif

#endif /* OPENSHIELDHIT_VERSION_H */
