#ifndef OPENSHIELDHIT_READLINE_H
#define OPENSHIELDHIT_READLINE_H

#include "openshieldhit/file.h"

#ifdef __cplusplus
extern "C" {
#endif

#define OSH_READLINE_COMMENT "#*"
#define OSH_MAX_LINE_LENGTH 4096

/**
 * @brief Read the next non-empty, non-comment line from an input file.
 *
 * @details
 * Leading/trailing whitespace and inline comments are stripped. The caller
 * owns `*line` and must free it when done.
 *
 * @returns Parsed line length, or -1 on EOF/error.
 */
int osh_readline(struct oshfile *oshf, char **line, int *lineno);

/**
 * @brief Split a mutable string into whitespace-delimited tokens in place.
 */
int osh_tokenise(char *line, char **words, int max_words);

/**
 * @brief Read the next non-comment line and split it into key plus arguments.
 *
 * @details
 * `*lline` is allocated by the function and must be freed by the caller.
 * `*kkey` and `*aargs` point into that same buffer.
 *
 * @returns Line length, or -1 on EOF/error.
 */
int osh_readline_key(struct oshfile *oshf, char **lline, char **kkey, char **aargs, int *lineno);

/**
 * @brief Lower-case a NUL-terminated string in-place.
 */
void osh_lower_inplace(char *s);

#ifdef __cplusplus
}
#endif

#endif /* OPENSHIELDHIT_READLINE_H */
