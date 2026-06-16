#include "osh_file.h"

#if defined(_WIN32)
#include <direct.h>
#include <errno.h>
#include <sys/stat.h>
#include <windows.h>
#define mkdir(path, mode) _mkdir(path)
#else
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#endif
#include <stdlib.h>
#include <string.h>

#include "openshieldhit/status.h"
#include "osh_abort.h"

static enum osh_status _mapfile(struct oshfile *oshf);
static enum osh_status _rewind_file(struct oshfile *oshf);
static int _is_sep(char c);
static int _mode_is_dir(int mode);
#if defined(_WIN32)
static int _is_drive_letter(char c);
#endif
static int _path_is_absolute(char const *path);
static size_t _path_root_len(char const *path);

struct oshfile *osh_fopen(char const *filename) {
    FILE *fp;
    struct oshfile *oshf;

    fp = fopen(filename, "r");
    if (!fp) {
        return NULL;
    }

    oshf = malloc(sizeof(struct oshfile));
    if (!oshf) {
        osh_abort_oomf("osh_fopen");
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
        osh_abort_oomf("osh_mapfile: failed to allocate memory for line map");
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

/**
 * @brief Create @p path and any missing parent directories.
 *
 * @details
 * The helper treats an already existing directory as success and rejects
 * existing non-directory paths. Absolute Windows drive roots such as `C:/`
 * and POSIX roots such as `/` are preserved and not created as intermediate
 * components.
 *
 * @param[in] path  Directory path to create.
 *
 * @returns OSH_OK on success, OSH_EINVAL for invalid input, OSH_ENOMEM for
 *          allocation failure, or OSH_EIO when filesystem operations fail.
 */
enum osh_status osh_path_ensure_dir(char const *path) {
    char *tmp;
    char *p;
    size_t len;
    size_t root_len;
    struct stat st;

    if (!path || !path[0]) {
        return OSH_EINVAL;
    }

    if (stat(path, &st) == 0) {
        return _mode_is_dir(st.st_mode) ? OSH_OK : OSH_EIO;
    }

    len = strlen(path);
    tmp = (char *) malloc(len + 1u);
    if (!tmp) {
        return OSH_ENOMEM;
    }
    memcpy(tmp, path, len + 1u);

    root_len = _path_root_len(tmp);
    for (p = tmp + root_len; *p; ++p) {
        if (!_is_sep(*p)) {
            continue;
        }
        if ((p > tmp) && _is_sep(*(p - 1))) {
            continue;
        }
        *p = '\0';
        if (tmp[0] != '\0' && stat(tmp, &st) != 0) {
            if (mkdir(tmp, 0777) != 0 && errno != EEXIST) {
                free(tmp);
                return OSH_EIO;
            }
        } else if (!_mode_is_dir(st.st_mode)) {
            free(tmp);
            return OSH_EIO;
        }
        *p = '/';
    }

    if (stat(tmp, &st) != 0) {
        if (mkdir(tmp, 0777) != 0 && errno != EEXIST) {
            free(tmp);
            return OSH_EIO;
        }
    } else if (!_mode_is_dir(st.st_mode)) {
        free(tmp);
        return OSH_EIO;
    }

    free(tmp);
    return OSH_OK;
}

enum osh_status osh_dir_foreach_file(char const *dir, osh_dir_iter_fn fn, void *user) {
    if (!dir || !fn) {
        return OSH_EINVAL;
    }

#if defined(_WIN32)
    {
        WIN32_FIND_DATAA fd;
        HANDLE hfind;
        char pattern[4096];
        char path[4096];

        if (snprintf(pattern, sizeof(pattern), "%s\\*", dir) < 0) {
            return OSH_EIO;
        }
        hfind = FindFirstFileA(pattern, &fd);
        if (hfind == INVALID_HANDLE_VALUE) {
            return OSH_EIO;
        }

        do {
            if ((fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0) {
                continue;
            }
            if (snprintf(path, sizeof(path), "%s\\%s", dir, fd.cFileName) < 0) {
                FindClose(hfind);
                return OSH_EIO;
            }
            if (!fn(path, user)) {
                break;
            }
        } while (FindNextFileA(hfind, &fd));

        FindClose(hfind);
        return OSH_OK;
    }
#else
    {
        DIR *d;
        struct dirent *ent;
        char path[4096];

        d = opendir(dir);
        if (!d) {
            return OSH_EIO;
        }

        while ((ent = readdir(d)) != NULL) {
            struct stat st;

            if (strcmp(ent->d_name, ".") == 0 || strcmp(ent->d_name, "..") == 0) {
                continue;
            }
            if (snprintf(path, sizeof(path), "%s/%s", dir, ent->d_name) < 0) {
                closedir(d);
                return OSH_EIO;
            }
            if (stat(path, &st) != 0) {
                closedir(d);
                return OSH_EIO;
            }
            if (!S_ISREG(st.st_mode)) {
                continue;
            }
            if (!fn(path, user)) {
                break;
            }
        }

        closedir(d);
        return OSH_OK;
    }
#endif
}

/**
 * @brief Return non-zero when @p mode describes a directory.
 *
 * @details
 * MSVC does not provide POSIX `S_ISDIR` consistently as a macro in all build
 * environments, so keep the check local and portable.
 */
static int _mode_is_dir(int mode) {
#if defined(_WIN32)
    return (mode & _S_IFMT) == _S_IFDIR;
#else
    return S_ISDIR(mode);
#endif
}

/**
 * @brief Return the rooted prefix length of @p path.
 *
 * @details
 * Returns the number of bytes that belong to the non-creatable root portion of
 * an absolute path. Examples: `/` -> 1, `C:/` -> 3, `//server/share/` -> the
 * length up to and including the share separator. Relative paths return 0.
 *
 * @param[in] path  Path string to inspect.
 *
 * @returns Root prefix length in bytes, or 0 for relative paths.
 */
static size_t _path_root_len(char const *path) {
#if defined(_WIN32)
    size_t i;
    int sep_count;
#endif

    if (!path || !path[0]) {
        return 0u;
    }

#if defined(_WIN32)
    if (_is_drive_letter(path[0]) && path[1] == ':' && _is_sep(path[2])) {
        return 3u;
    }
    if (_is_sep(path[0]) && _is_sep(path[1])) {
        i = 2u;
        sep_count = 0;
        while (path[i]) {
            if (_is_sep(path[i])) {
                sep_count++;
                if (sep_count == 2) {
                    return i + 1u;
                }
            }
            i++;
        }
        return 2u;
    }
#endif

    if (_is_sep(path[0])) {
        return 1u;
    }

    return 0u;
}
