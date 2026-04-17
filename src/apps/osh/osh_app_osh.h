#ifndef OSH_APP_OSH_H
#define OSH_APP_OSH_H

#include "material/osh_material.h"
#include "openshieldhit/beam.h"
#include "openshieldhit/logger.h"
#include "openshieldhit/status.h"

int osh_beam_setup_from_path(char const *path, struct osh_logger *lg, struct osh_beam_workspace **wb_out);
enum osh_status
osh_material_setup_from_path(char const *path, struct osh_logger *lg, struct material_workspace **wm_out);

#endif /* OSH_APP_OSH_H */
