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

struct oshfile *osh_fopen(const char *filename) {
    FILE *fp;
    struct oshfile *oshf;

    fp = fopen(filename, "r");
    if (!fp) {
        osh_error(EX_IOERR, "Could not open file: %s", filename);
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

    _mapfile(oshf);
    if (oshf->map_len < 1) {
        osh_error(EX_IOERR, "osh_fopen: file appears to contain no newlines.");
    }

    return oshf;
}

void osh_fclose(struct oshfile *oshf) {
    if (!oshf)
        return;

    if (oshf->fp)
        fclose(oshf->fp);
    free(oshf->filename);
    free(oshf->map);
    free(oshf);
}

int osh_file_lineno(const struct oshfile *oshf) {
    long int pos;
    int low = 0, high, mid;

    if (!oshf || !oshf->fp || !oshf->map || oshf->map_len == 0)
        return -1;

    pos = ftell(oshf->fp);
    if (pos < 0)
        return -1;

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
    char c;
    long int i;

    if (!oshf || !oshf->fp) {
        osh_error(EX_SOFTWARE, "osh_mapfile: null file pointer");
    }

    /* Rewind and count newlines */
    rewind(oshf->fp);
    oshf->map_len = 0;

    for (c = getc(oshf->fp); c != EOF; c = getc(oshf->fp)) {
        if (c == '\n') {
            oshf->map_len++;
        }
    }

    if (oshf->map_len < 1) {
        osh_warn("osh_mapfile: file appears to contain no newlines.");
    }

    /* Allocate and fill the map with byte address of each new line */
    oshf->map = calloc(oshf->map_len, sizeof(long int));
    if (!oshf->map) {
        osh_alloc_failed("osh_mapfile: failed to allocate memory for line map");
    }

    rewind(oshf->fp);

    i = 0;
    for (c = getc(oshf->fp); c != EOF; c = getc(oshf->fp)) {
        if (c == '\n') {
            oshf->map[i] = ftell(oshf->fp);
            i++;
        }
    }

    rewind(oshf->fp);

    return 1; /* Success */
}
