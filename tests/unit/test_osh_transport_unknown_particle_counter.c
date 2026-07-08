#include <string.h>

#include "common/osh_particle_pool.h"
#include "common/osh_step_segment.h"
#include "material/runtime/osh_material_runtime.h"
#include "particle/osh_particle.h"
#include "test_assert.h"
#include "transport/osh_transport.h"
#include "transport/osh_transport_ion_step.h"

static void test_unknown_particle_drop_is_counted(void) {
    struct osh_particle_pool pool;
    struct osh_transport_context transport_ctx;
    struct osh_material_runtime material_rt;
    struct osh_step_segment segment;
    struct particle unknown;
    unsigned int projectile_z[1];
    enum osh_status rc;

    memset(&transport_ctx, 0, sizeof(transport_ctx));
    memset(&material_rt, 0, sizeof(material_rt));
    memset(&segment, 0, sizeof(segment));
    memset(&unknown, 0, sizeof(unknown));

    projectile_z[0] = 1u;
    material_rt.projectile_z = projectile_z;
    material_rt.nprojectiles = 1u;

    unknown.mass = 3727.0;
    unknown.charge = 2;
    unknown.z = 2u;
    unknown.a = 4u;
    unknown.is_nucleus = 1u;

    ASSERT_TRUE(osh_particle_pool_init(&pool, 1u) == OSH_OK);
    pool.n = 1u;
    pool.e[0] = 10.0;
    pool.ux[0] = 1.0;
    pool.species[0] = &unknown;

    rc = osh_transport_ion_step(
        &pool, 0u, NULL, &segment, 1u, NULL, &transport_ctx, NULL, &material_rt, NULL, NULL, NULL);

    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(pool.e[0] == 0.0);
    ASSERT_TRUE(pool.n_unknown_dropped == 1u);
    ASSERT_TRUE(pool.n_dropped == 0u);

    osh_particle_pool_free(&pool);
}

int main(void) {
    test_unknown_particle_drop_is_counted();
    return 0;
}
