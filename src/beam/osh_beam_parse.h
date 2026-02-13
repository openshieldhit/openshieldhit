#ifndef OSH_BEAM_PARSE_H
#define OSH_BEAM_PARSE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "beam/osh_beam.h"
#include "common/osh_file.h"

int osh_beam_parse(struct oshfile *oshf, struct beam_workspace *beam);

#endif /* OSH_BEAM_PARSE_H */