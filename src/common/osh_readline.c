#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "osh_logger.h"
#include "osh_readline.h"

static int _is_comment(char c);

int osh_readline(struct oshfile *oshf, char **line, int *lineno) {
    char buff[OSH_MAX_LINE_LENGTH];
    int pos0, pos1;
    int slen;
    int i;

    if (!oshf || !oshf->fp) {
        osh_error(EX_SOFTWARE, "osh_readline: file pointer is NULL");
    }

    if (*line == NULL) {
        *line = calloc(OSH_MAX_LINE_LENGTH, sizeof(char));
        if (!*line) {
            osh_alloc_failed(OSH_MAX_LINE_LENGTH * sizeof(char));
        }
    }

    while (fgets(buff, OSH_MAX_LINE_LENGTH, oshf->fp)) {
        oshf->lineno++;
        slen = strlen(buff);
        pos0 = 0;
        pos1 = -1;

        for (i = 0; i < slen; ++i) {
            if (_is_comment(buff[i]) || buff[i] == '\0') {
                break;
            }
            if (!isspace((unsigned char)buff[i])) {
                if (pos1 == -1)
                    pos0 = i;
                pos1 = i;
            }
        }

        if (pos1 >= pos0) {
            int len = pos1 - pos0 + 1;
            snprintf(*line, len + 1, "%s", buff + pos0);
            if (lineno)
                *lineno = oshf->lineno;
            return len;
        }
    }

    return -1; /* EOF or error */
}

int osh_readline_key(struct oshfile *oshf, char **lline, char **kkey, char **aargs, int *lineno) {
    char buff[OSH_MAX_LINE_LENGTH];
    char *line = NULL;
    int slen;
    int key_end;
    int i;
    int j;
    int args_start;
    int last_char;

    while (fgets(buff, OSH_MAX_LINE_LENGTH, oshf->fp)) {
        oshf->lineno++;
        slen = strlen(buff);

        /* Skip empty lines and comments */
        i = 0;
        while (isspace((unsigned char)buff[i]))
            i++;
        if (buff[i] == '\0' || _is_comment(buff[i]))
            continue;

        /* Allocate and copy line */
        line = calloc(OSH_MAX_LINE_LENGTH, sizeof(char));
        if (!line) {
            osh_alloc_failed(OSH_MAX_LINE_LENGTH * sizeof(char));
            return -1;
        }
        strncpy(line, buff, OSH_MAX_LINE_LENGTH);

        /* Find key */
        i = 0;
        while (isspace((unsigned char)line[i]))
            i++;
        *kkey = line + i;
        key_end = i;
        while (line[key_end] && !isspace((unsigned char)line[key_end]) && !_is_comment(line[key_end]))
            key_end++;
        if (line[key_end])
            line[key_end] = '\0';

        /* Find args */
        args_start = key_end + 1;
        while (line[args_start] && isspace((unsigned char)line[args_start]))
            args_start++;
        if (line[args_start] == '\0' || _is_comment(line[args_start])) {
            *aargs = NULL;
        } else {
            *aargs = line + args_start;
            j = args_start;
            last_char = j;
            while (line[j] && !_is_comment(line[j])) {
                if (!isspace((unsigned char)line[j]))
                    last_char = j;
                j++;
            }
            line[last_char + 1] = '\0';
        }

        *lline = line;
        if (lineno)
            *lineno = oshf->lineno;
        return slen;
    }

    *lline = NULL;
    *kkey = NULL;
    *aargs = NULL;
    return -1;
}

/**
 * @brief Checks if character c is a comment marker
 *
 * @param[in]  c - character to be checked
 *
 * @returns 1 if it is a comment, 0 if it is not a comment
 *
 * @author Niels Bassler
 */
static int _is_comment(char c) {

    int i = 0;
    int cl;

    cl = strlen(OSH_READLINE_COMMENT);

    for (i = 0; i < cl; i++) {
        if (c == OSH_READLINE_COMMENT[i])
            return 1;
    }
    return 0;
}
