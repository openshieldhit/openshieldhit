#ifndef OPENSHIELDHIT_SIMULATION_H
#define OPENSHIELDHIT_SIMULATION_H

#include "openshieldhit/beam.h"
#include "openshieldhit/geometry.h"
#include "openshieldhit/material.h"
#include "openshieldhit/scoring.h"
#include "openshieldhit/status.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Opaque handle for a compiled, ready-to-run simulation.
 *
 * @details
 * Holds all four runtime representations compiled from the cold workspaces.
 * The cold workspaces are borrowed (not owned); the caller must keep them
 * alive for the lifetime of the simulation object.
 */
struct osh_simulation;

/**
 * @brief Compile four cold workspaces into a simulation ready to run.
 *
 * @details
 * Performs zone-to-material index resolution, compiles geometry and scoring
 * runtimes, prepares transport tables, and initialises the beam source.
 *
 * @param[in]  beam     Prepared beam workspace.
 * @param[in]  geo      Prepared geometry workspace.
 * @param[in]  mat      Finalized material workspace.
 * @param[in]  scoring  Parsed scoring workspace.
 * @param[out] sim_out  Receives the new simulation handle on success.
 *
 * @returns OSH_OK on success, or an error code.
 */
enum osh_status osh_simulation_create(struct osh_beam_workspace *beam,
                                      struct osh_geometry_workspace *geo,
                                      struct osh_material_workspace *mat,
                                      struct osh_scoring_workspace *scoring,
                                      struct osh_simulation **sim_out);

/**
 * @brief Run the simulation and save outputs to @p out_dir.
 *
 * @details
 * Drives the transport loop, postprocesses scoring accumulators, and writes
 * all output files.
 *
 * @param[in] sim      Simulation handle created by osh_simulation_create().
 * @param[in] out_dir  Directory for output files.
 *
 * @returns OSH_OK on success, or an error code.
 */
enum osh_status osh_simulation_run(struct osh_simulation *sim, char const *out_dir);

/**
 * @brief Release the simulation.
 *
 * @details
 * Frees all runtime resources owned by the simulation.  The cold workspaces
 * passed to osh_simulation_create() are not freed.
 *
 * @param[in] sim  May be NULL (no-op).
 *
 * @returns OSH_OK always.
 */
enum osh_status osh_simulation_free(struct osh_simulation *sim);

#ifdef __cplusplus
}
#endif

#endif /* OPENSHIELDHIT_SIMULATION_H */
