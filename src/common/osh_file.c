#include "osh_file.h"

#include <stdlib.h>
#include <string.h>

#include "osh_logger.h"

/**
 * @brief Creates a byte-map for a given file *oshf.
 *
 * @details This function reads the file line by line, counts the number of
 * newlines, and creates a map of byte offsets for each line in the file. The
 * file pointer is rewound after this function is called.
 *
 * @param[in,out] oshf Pointer to struct oshfile struct.
 *
 * @returns 1 on success, or exits with an error if the file pointer is NULL or
 * memory allocation fails.
 *
 * @author Niels Bassler
 */
static int _mapfile(struct oshfile *oshf);
static int _rewind_file(struct oshfile *oshf);

struct oshfile *osh_fopen(char const *filename) {
    FILE *fp;
    struct oshfile *oshf;

    fp = fopen(filename, "r");
    if (!fp) {
        osh_error("Could not open file: %s", filename);
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

    if (!_mapfile(oshf)) {
        osh_fclose(oshf);
        return NULL;
    }
    if (oshf->map_len < 1) {
        osh_error("osh_fopen: file appears to contain no newlines.");
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

static int _mapfile(struct oshfile *oshf) {
    int c;
    long int i;
    size_t map_slots;

    if (!oshf || !oshf->fp) {
        osh_error("osh_mapfile: null file pointer");
        return 0;
    }

    /* Rewind and count newlines */
    if (!_rewind_file(oshf)) {
        return 0;
    }
    oshf->map_len = 0;

    for (c = fgetc(oshf->fp); c != EOF; c = fgetc(oshf->fp)) {
        if (c == '\n') {
            oshf->map_len++;
        }
    }
    if (ferror(oshf->fp)) {
        osh_error("osh_mapfile: failed while scanning file '%s'", oshf->filename);
        return 0;
    }

    if (oshf->map_len < 1) {
        osh_warn("osh_mapfile: file appears to contain no newlines.");
    }

    /* Allocate and fill the map with byte address of each new line */
    map_slots = (oshf->map_len > 0) ? (size_t) oshf->map_len : 1u;
    oshf->map = calloc(map_slots, sizeof(long int));
    if (!oshf->map) {
        osh_alloc_failed("osh_mapfile: failed to allocate memory for line map");
    }

    if (!_rewind_file(oshf)) {
        return 0;
    }

    i = 0;
    for (c = fgetc(oshf->fp); c != EOF; c = fgetc(oshf->fp)) {
        if (c == '\n') {
            oshf->map[i] = ftell(oshf->fp);
            i++;
        }
    }
    if (ferror(oshf->fp)) {
        osh_error("osh_mapfile: failed while building line map for '%s'", oshf->filename);
        return 0;
    }

    if (!_rewind_file(oshf)) {
        return 0;
    }

    return 1; /* Success */
}

static int _rewind_file(struct oshfile *oshf) {
    if (fseek(oshf->fp, 0L, SEEK_SET) != 0) {
        osh_error("osh_mapfile: failed to rewind file '%s'", oshf->filename);
        return 0;
    }
    clearerr(oshf->fp);
    oshf->lineno = 0;
    return 1;
}

int osh_relative_path_to_file(char **out, char const *base_dir, char const *rel_path) {
    size_t blen;
    size_t rlen;
    char *result;

    if (!out || !rel_path) {
        return -1;
    }

    if (rel_path[0] == '/') {
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
    char const *sep;
    size_t len;
    char *dir;

    if (!path) {
        return NULL;
    }
    sep = strrchr(path, '/');
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
