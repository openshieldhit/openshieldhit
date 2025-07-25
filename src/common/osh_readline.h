#ifndef OSH_READLINE_H
#define OSH_READLINE_H

#include "osh_file.h"

#define OSH_READLINE_COMMENT "#*!"  /* list of characters which will be treated as comment prefix */
#define OSH_MAX_LINE_LENGTH 4096    /* maximum line length inside a file to be read in bytes */

/**
 * @brief Reads the next non-empty, non-comment line from a osh_file.
 *
 * @details Strips leading/trailing whitespace and inline comments.
 * Caller must free *line after use.
 *
 * @param[in,out] oshf  Open file handle
 * @param[out]    line  Pointer to allocated line buffer (allocated by function)
 * @param[out]    lineno Optional pointer to receive line number (can be NULL)
 * @returns Length of parsed line (excluding null terminator), or -1 on EOF
 *
 * @author Niels Bassler
 */
int osh_readline(struct oshfile *oshf, char **line, int *lineno);

/**
 * @brief Reads the next non-comment line, splits into key and arguments.
 *
 * Allocates memory for the line; caller must free *lline unless -1 is returned.
 * *kkey and *aargs point within *lline; *aargs is NULL if no arguments.
 *
 * @param[in]  oshf    Input file handle.
 * @param[out] lline   Pointer to allocated line buffer.
 * @param[out] kkey    Pointer to key (first word).
 * @param[out] aargs   Pointer to arguments, or NULL.
 * @param[out] lineno  Line number of returned line.
 *
 * @return Length of line (excluding null byte), or -1 on error/EOF.
 */
int osh_readline_key(struct oshfile *oshf, char **lline, char **kkey, char **aargs, int *lineno);

#endif /* OSH_READLINE_H */