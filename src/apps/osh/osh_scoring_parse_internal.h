#ifndef OSH_SCORING_PARSE_INTERNAL_H
#define OSH_SCORING_PARSE_INTERNAL_H

#include "scoring/osh_scoring.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Section-level line dispatchers.
 *
 * Each function receives one already-tokenised line (words[0] is the
 * lower-cased keyword) and updates the section being accumulated.
 * @p found_out is set to 1 when the keyword was recognised, 0 otherwise.
 * An unrecognised keyword is not an error — the caller decides whether to warn.
 */

/**
 * @brief Dispatch one tokenised line into a Filter section.
 */
enum osh_status osh_scoring_parse_filter_line(struct osh_scoring_filter_def *fil,
                                              char **words,
                                              int nwords,
                                              char const *path,
                                              unsigned int lineno,
                                              int *found_out);

/**
 * @brief Dispatch one tokenised line into a Settings section.
 */
enum osh_status osh_scoring_parse_settings_line(struct osh_scoring_settings_def *set,
                                                char **words,
                                                int nwords,
                                                char const *path,
                                                unsigned int lineno,
                                                int *found_out);

/**
 * @brief Dispatch one tokenised line into a Geometry section.
 */
enum osh_status osh_scoring_parse_geometry_line(struct osh_scoring_geometry_def *geo,
                                                char **words,
                                                int nwords,
                                                char const *path,
                                                unsigned int lineno,
                                                int *found_out);

/**
 * @brief Dispatch one tokenised line into an Output section.
 */
enum osh_status osh_scoring_parse_output_line(struct osh_scoring_output_def *out,
                                              char **words,
                                              int nwords,
                                              char const *path,
                                              unsigned int lineno,
                                              int *found_out);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SCORING_PARSE_INTERNAL_H */
