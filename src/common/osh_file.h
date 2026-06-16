#ifndef OSH_FILE_H
#define OSH_FILE_H

#include <stddef.h>
#include <stdio.h>

#include "openshieldhit/status.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Lightweight text-file wrapper with filename and line metadata.
 *
 * @details
 * Parser-oriented helper used by app-layer importers that need stable
 * filename/line-number diagnostics while reading ASCII input files.
 */
struct oshfile {
    FILE *fp;
    long int *map;
    char *filename;
    size_t map_size;
    int map_len;
    int lineno;
};

struct oshfile *osh_fopen(char const *filename);
void osh_fclose(struct oshfile *oshf);
int osh_file_lineno(struct oshfile const *oshf);

/**
 * @brief Resolve a referenced path relative to the containing input file.
 *
 * @details
 * If @p rel_path is already absolute it is copied as-is. Otherwise it is
 * resolved relative to @p base_dir. On Windows builds, drive-absolute paths
 * such as `C:/...` and `C:\\...` are recognized as absolute; bare
 * drive-relative paths such as `C:foo` are not.
 *
 * The returned string is heap-allocated and owned by the caller.
 *
 * @returns 0 on success, -1 on allocation or argument failure.
 */
int osh_relative_path_to_file(char **out, char const *base_dir, char const *rel_path);

/**
 * @brief Normalize path separators to '/' in-place.
 *
 * @details
 * On Windows this converts `\\` to `/`. On non-Windows platforms it is a
 * no-op. Safe to call with NULL.
 */
void osh_path_normalize(char *path);

/**
 * @brief Return the directory portion of a path as a new string.
 *
 * @details
 * Accepts both '/' and '\\' as separators so callers remain robust when a
 * Windows-style path reaches this helper before normalization.
 *
 * @returns Newly allocated directory string, or NULL on missing separator or
 * allocation failure.
 */
char *osh_path_dirname(char const *path);

/**
 * @brief Create a directory path and any missing parent directories.
 *
 * @details
 * Existing directories are treated as success. The helper accepts both `/`
 * and `\\` separators and handles rooted paths on the current platform.
 *
 * @return OSH_OK on success, OSH_EINVAL for invalid arguments, OSH_ENOMEM for
 *         allocation failure, or OSH_EIO on filesystem errors.
 */
enum osh_status osh_path_ensure_dir(char const *path);

/**
 * @brief Callback used by osh_dir_foreach_file().
 *
 * @param[in] path  Full path to one non-directory entry.
 * @param[in] user  Caller context.
 *
 * @return Non-zero to continue, zero to stop iteration early.
 */
typedef int (*osh_dir_iter_fn)(char const *path, void *user);

/**
 * @brief Iterate over non-directory entries in a directory.
 *
 * @details
 * This helper centralizes the platform-specific directory walking needed by
 * callers such as the DICOM CT series loader. The callback receives one full
 * path per non-directory entry and may stop iteration early by returning 0.
 *
 * @return OSH_OK on success, OSH_EIO if the directory cannot be opened or
 * walked, or OSH_EINVAL for invalid arguments.
 */
enum osh_status osh_dir_foreach_file(char const *dir, osh_dir_iter_fn fn, void *user);

#ifdef __cplusplus
}
#endif

#endif /* OSH_FILE_H */
