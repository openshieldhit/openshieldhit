#ifndef OPENSHIELDHIT_H
#define OPENSHIELDHIT_H

/*
 * Convenience umbrella header for embedders who want the most common public
 * OpenShieldHIT interfaces in one include.
 *
 * Header roles:
 *   openshieldhit/openshieldhit.h
 *     High-level opaque context API for whole-run orchestration.
 *
 *   openshieldhit/beam.h
 *     Public cold beam data model owned by `osh_core`.
 *
 * Application / parser code lives under src/apps/... and is intentionally
 * not part of the installed public API.
 */

#include "openshieldhit/beam.h"
#include "openshieldhit/beam_defs.h"
#include "openshieldhit/file.h"
#include "openshieldhit/logger.h"
#include "openshieldhit/openshieldhit.h"
#include "openshieldhit/readline.h"
#include "openshieldhit/status.h"

#endif /* OPENSHIELDHIT_H */
