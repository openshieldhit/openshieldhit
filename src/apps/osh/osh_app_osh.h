#ifndef OSH_APP_OSH_H
#define OSH_APP_OSH_H

#include "openshieldhit/beam.h"
#include "openshieldhit/geometry.h"
#include "openshieldhit/logger.h"
#include "openshieldhit/material.h"
#include "openshieldhit/scoring.h"
#include "openshieldhit/status.h"

enum osh_status osh_beam_setup_from_path(char const *path, struct osh_logger *lg, struct osh_beam_workspace **wb_out);
enum osh_status
osh_geometry_setup_from_path(char const *path, struct osh_logger *lg, struct osh_geometry_workspace **ws_out);
enum osh_status
osh_material_setup_from_path(char const *path, struct osh_logger *lg, struct osh_material_workspace **wm_out);
enum osh_status
osh_scoring_setup_from_path(char const *path, struct osh_logger *lg, struct osh_scoring_workspace **ws_out);

#endif /* OSH_APP_OSH_H */
