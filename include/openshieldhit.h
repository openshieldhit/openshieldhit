#ifndef OPENSHIELDHIT_H
#define OPENSHIELDHIT_H

/**
 * @file openshieldhit.h
 * @brief Convenience umbrella header — includes all public OpenShieldHIT interfaces.
 *
 * @details
 * Including this header is equivalent to including each public header
 * individually.  Prefer explicit includes in library and adapter code so
 * that dependencies are visible; use this umbrella in application code
 * where conciseness matters more.
 *
 * Application / parser code lives under src/apps/... and is intentionally
 * not part of the installed public API.
 */

#include "openshieldhit/beam.h"
#include "openshieldhit/beam_defs.h"
#include "openshieldhit/file.h"
#include "openshieldhit/geometry.h"
#include "openshieldhit/geometry_defs.h"
#include "openshieldhit/logger.h"
#include "openshieldhit/material.h"
#include "openshieldhit/readline.h"
#include "openshieldhit/scoring.h"
#include "openshieldhit/status.h"
#include "openshieldhit/version.h"

#endif /* OPENSHIELDHIT_H */
