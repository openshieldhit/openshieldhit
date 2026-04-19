#include "osh_file.h"

#include <stdlib.h>
#include <string.h>

#include "openshieldhit/status.h"
#include "osh_logger.h"

static enum osh_status _mapfile(struct oshfile *oshf);
static enum osh_status _rewind_file(struct oshfile *oshf);
static int _is_sep(char c);
#if defined(_WIN32)
static int _is_drive_letter(char c);
#endif
static int _path_is_absolute(char const *path);

struct oshfile *osh_fopen(char const *filename) {
    FILE *fp;
    struct oshfile *oshf;

    fp = fopen(filename, "r");
    if (!fp) {
        return NULL;
    }

    oshf = malloc(sizeof(struct oshfile));
    if (!oshf) {
        osh_alloc_failed("osh_fopen");
    }

    oshf->fp = fp;
    oshf->filename = strdup(filename);
    oshf->lineno = 0; /* current line number */

    oshf->map = NULL;
    oshf->map_len = 0; /* number of lines (entries) of the map */

    if (_mapfile(oshf) != OSH_OK) {
        osh_fclose(oshf);
        return NULL;
    }
    if (oshf->map_len < 1) {
        osh_fclose(oshf);
        return NULL;
    }

    return oshf;
}

void osh_fclose(struct oshfile *oshf) {
    if (!oshf) {
        return;
    }

    if (oshf->fp) {
        fclose(oshf->fp);
    }
    free(oshf->filename);
    free(oshf->map);
    free(oshf);
}

int osh_file_lineno(const struct oshfile *oshf) {
    long int pos;
    int low = 0;
    int high;
    int mid;

    if (!oshf || !oshf->fp || !oshf->map || oshf->map_len == 0) {
        return -1;
    }

    pos = ftell(oshf->fp);
    if (pos < 0) {
        return -1;
    }

    high = oshf->map_len - 1;

    while (low <= high) {
        mid = (low + high) / 2;
        if (pos <= oshf->map[mid]) {
            high = mid - 1;
        } else {
            low = mid + 1;
        }
    }

    return low + 1; /* Line numbers are 1-based */
}

/**
 * @brief Build the newline offset map used for line-number queries.
 *
 * @details
 * Scans the file twice: first to count newline characters, then to record the
 * byte offset immediately after each newline. The file is rewound before and
 * after the scan so callers can continue reading from the beginning.
 *
 * @param[in,out] oshf  File wrapper whose `map` and `map_len` fields are filled.
 *
 * @returns OSH_OK on success, or an OSH_E* code on failure.
 */
static enum osh_status _mapfile(struct oshfile *oshf) {
    int c;
    long int i;
    size_t map_slots;

    if (!oshf || !oshf->fp) {
        return OSH_EINVAL;
    }

    /* Rewind and count newlines */
    if (_rewind_file(oshf) != OSH_OK) {
        return OSH_EIO;
    }
    oshf->map_len = 0;

    for (c = fgetc(oshf->fp); c != EOF; c = fgetc(oshf->fp)) {
        if (c == '\n') {
            oshf->map_len++;
        }
    }
    if (ferror(oshf->fp)) {
        return OSH_EIO;
    }

    /* Allocate and fill the map with byte address of each new line */
    map_slots = (oshf->map_len > 0) ? (size_t) oshf->map_len : 1u;
    oshf->map = calloc(map_slots, sizeof(long int));
    if (!oshf->map) {
        osh_alloc_failed("osh_mapfile: failed to allocate memory for line map");
    }

    if (_rewind_file(oshf) != OSH_OK) {
        return OSH_EIO;
    }

    i = 0;
    for (c = fgetc(oshf->fp); c != EOF; c = fgetc(oshf->fp)) {
        if (c == '\n') {
            oshf->map[i] = ftell(oshf->fp);
            i++;
        }
    }
    if (ferror(oshf->fp)) {
        return OSH_EIO;
    }

    if (_rewind_file(oshf) != OSH_OK) {
        return OSH_EIO;
    }

    return OSH_OK;
}

static enum osh_status _rewind_file(struct oshfile *oshf) {
    if (!oshf || !oshf->fp) {
        return OSH_EINVAL;
    }
    if (fseek(oshf->fp, 0L, SEEK_SET) != 0) {
        return OSH_EIO;
    }
    clearerr(oshf->fp);
    oshf->lineno = 0;
    return OSH_OK;
}

/**
 * @brief Return whether @p c is a recognized path separator.
 *
 * @details
 * Accept both '/' and '\\' so helpers remain robust even when a caller forgot
 * to normalize Windows-style paths first.
 */
static int _is_sep(char c) {
    return (c == '/') || (c == '\\');
}

#if defined(_WIN32)
/**
 * @brief Return whether @p c is an ASCII drive letter.
 */
static int _is_drive_letter(char c) {
    return ((c >= 'A') && (c <= 'Z')) || ((c >= 'a') && (c <= 'z'));
}
#endif

/**
 * @brief Return whether @p path is already absolute for the current platform.
 *
 * @details
 * On all platforms, a leading separator denotes an absolute/rooted path.
 * On Windows, also recognize drive-absolute paths (`C:/...`, `C:\\...`) and
 * UNC paths (`\\\\server\\share`, `//server/share`). A bare drive-relative
 * path like `C:foo` is intentionally not treated as absolute.
 */
static int _path_is_absolute(char const *path) {
    if (!path || !path[0]) {
        return 0;
    }

    if (_is_sep(path[0])) {
        return 1;
    }

#if defined(_WIN32)
    if (_is_drive_letter(path[0]) && path[1] == ':' && _is_sep(path[2])) {
        return 1;
    }
#endif

    return 0;
}

int osh_relative_path_to_file(char **out, char const *base_dir, char const *rel_path) {
    size_t blen;
    size_t rlen;
    char *result;

    if (!out || !rel_path) {
        return -1;
    }

    if (_path_is_absolute(rel_path)) {
        /* already absolute — copy as-is */
        rlen = strlen(rel_path);
        result = (char *) malloc(rlen + 1);
        if (!result) {
            return -1;
        }
        memcpy(result, rel_path, rlen + 1);
        *out = result;
        return 0;
    }

    blen = base_dir ? strlen(base_dir) : 0;
    rlen = strlen(rel_path);
    result = (char *) malloc(blen + 1 + rlen + 1); /* base + '/' + rel + NUL */
    if (!result) {
        return -1;
    }
    if (blen > 0) {
        memcpy(result, base_dir, blen);
        result[blen] = '/';
        memcpy(result + blen + 1, rel_path, rlen + 1);
    } else {
        memcpy(result, rel_path, rlen + 1);
    }
    *out = result;
    return 0;
}

void osh_path_normalize(char *path) {
#ifdef _WIN32
    if (!path) {
        return;
    }
    for (; *path; ++path) {
        if (*path == '\\') {
            *path = '/';
        }
    }
#else
    (void) path; /* no-op on non-Windows */
#endif
}

char *osh_path_dirname(char const *path) {
    char const *p;
    char const *sep;
    size_t len;
    char *dir;

    if (!path) {
        return NULL;
    }
    sep = NULL;
    for (p = path; *p; ++p) {
        if (_is_sep(*p)) {
            sep = p;
        }
    }
    if (!sep || sep == path) {
        return NULL;
    }
    len = (size_t) (sep - path);
    dir = (char *) malloc(len + 1);
    if (!dir) {
        return NULL;
    }
    memcpy(dir, path, len);
    dir[len] = '\0';
    return dir;
}
