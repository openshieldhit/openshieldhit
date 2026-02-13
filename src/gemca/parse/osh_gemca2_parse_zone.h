#ifndef _OSH_GEMCA2_PARSE_ZONE
#define _OSH_GEMCA2_PARSE_ZONE

#include "common/osh_readline.h"
#include "gemca/osh_gemca2.h"

int osh_gemca_zone_init(struct zone **zone);
size_t osh_gemca_parse_count_zones(struct oshfile *shf);
int osh_gemca_parse_zones(struct oshfile *shf, struct gemca_workspace *g);

#endif /* _OSH_GEMCA2_PARSE_ZONE */
