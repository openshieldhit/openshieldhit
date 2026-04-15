#include "transport/osh_transport.h"
#include "transport/osh_transport_ion.h"

/*
 * Particle-type dispatcher for osh_transport_run_minimal().
 *
 * The public API dispatches here; each particle class has its own transport
 * module:
 *
 *   ions    → osh_transport_ion.c   (CSDA + Highland MCS + Bohr straggling)
 *   neutrons → osh_transport_neutron.c  (stub)
 *   photons  → osh_transport_photon.c   (stub)
 *
 * The beam->particle field (or PDG code) will drive dispatch once multiple
 * particle types are supported.  For now all transport goes to the ion loop.
 */
enum osh_status osh_transport_run_minimal(struct beam_workspace const *beam,
                                          struct gemca_runtime const *geom_rt,
                                          struct material_workspace const *materials,
                                          struct osh_material_runtime const *tables,
                                          struct osh_scoring_runtime *scoring) {
    return osh_transport_ion_run_minimal(beam, geom_rt, materials, tables, scoring);
}
