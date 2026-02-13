#ifndef _OSH_GEMCA2_PARSE_BODY
#define _OSH_GEMCA2_PARSE_BODY

#include "common/osh_readline.h"
#include "gemca/osh_gemca2.h"

int osh_gemca_body_init(struct body **body);
size_t osh_gemca_parse_count_bodies(struct oshfile *shf);
int osh_gemca_parse_bodies(struct oshfile *shf, struct gemca_workspace *g);

#endif /* _OSH_GEMCA2_PARSE_BODY */
