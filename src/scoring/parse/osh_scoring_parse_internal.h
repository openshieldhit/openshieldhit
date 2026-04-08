#ifndef OSH_SCORING_PARSE_INTERNAL_H
#define OSH_SCORING_PARSE_INTERNAL_H

#include "scoring/osh_scoring.h"

#ifdef __cplusplus
extern "C" {
#endif

enum osh_status osh_scoring_parse_filter_line(struct osh_scoring_filter_def *fil,
                                              char **words,
                                              int nwords,
                                              char const *path,
                                              unsigned int lineno,
                                              int *found_out);

enum osh_status osh_scoring_parse_settings_line(struct osh_scoring_settings_def *set,
                                                char **words,
                                                int nwords,
                                                char const *path,
                                                unsigned int lineno,
                                                int *found_out);

enum osh_status osh_scoring_parse_geometry_line(struct osh_scoring_geometry_def *geo,
                                                char **words,
                                                int nwords,
                                                char const *path,
                                                unsigned int lineno,
                                                int *found_out);

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
