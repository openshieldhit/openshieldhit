#include "transport/osh_transport_photon.h"

enum osh_status osh_transport_photon_run_minimal(struct beam_workspace const *beam,
                                                 struct gemca_runtime const *geom_rt,
                                                 struct material_workspace const *materials,
                                                 struct osh_material_runtime const *tables,
                                                 struct osh_scoring_runtime *scoring) {
    (void) beam;
    (void) geom_rt;
    (void) materials;
    (void) tables;
    (void) scoring;
    return OSH_ENOTSUP;
}
