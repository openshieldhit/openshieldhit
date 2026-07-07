#include <math.h>
#include <string.h>

#include "common/osh_particle_pool.h"
#include "common/osh_step_segment.h"
#include "material/runtime/osh_material_runtime.h"
#include "particle/osh_particle_pdg.h"
#include "random/osh_rng.h"
#include "scoring/runtime/osh_scoring_runtime.h"
#include "test_assert.h"
#include "transport/osh_transport.h"
#include "transport/osh_transport_ion.h"
#include "transport/osh_transport_ion_step.h"

/*
 * Argument-validation guards for the transport entry points.  These exercise the
 * early OSH_EINVAL returns without building a full runtime: the borrowed runtime
 * pointers are non-NULL sentinels that the guards reject on before ever
 * dereferencing them, so a valid, never-read address is all that is required.
 */

static void test_transport_run_rejects_bad_args(void) {
    struct osh_transport_context ctx;

    /* NULL context. */
    ASSERT_TRUE(osh_transport_run(NULL, NULL, NULL, NULL, NULL) == OSH_EINVAL);

    /* nstat == 0 is rejected before any pool is touched, so a zeroed context
     * (NULL pools) is enough to reach and return from that guard. */
    memset(&ctx, 0, sizeof(ctx));
    ctx.params.nstat = 0u;
    ASSERT_TRUE(osh_transport_run(&ctx, NULL, NULL, NULL, NULL) == OSH_EINVAL);
}

static void test_ion_run_range_validates_arguments(void) {
    struct osh_transport_context ctx;
    double storage = 0.0;
    void *nn = &storage; /* non-NULL sentinel; the guards never dereference it */
    size_t completed = 12345u;

    /* NULL context: rejected, and completed_out is zeroed up front. */
    ASSERT_TRUE(osh_transport_ion_run_range(NULL, nn, nn, nn, nn, NULL, 0u, 10u, &completed) == OSH_EINVAL);
    ASSERT_TRUE(completed == 0u);

    /* Missing pool / scratch on the context. */
    memset(&ctx, 0, sizeof(ctx));
    ASSERT_TRUE(osh_transport_ion_run_range(&ctx, nn, nn, nn, nn, NULL, 0u, 10u, &completed) == OSH_EINVAL);

    /* Wire non-NULL sentinel pool + scratch so the later guards are reachable. */
    ctx.ion_pool = nn;
    ctx.zone_refs = nn;
    ctx.dist_batch = nn;

    /* Empty / inverted range (hist_hi <= hist_lo). */
    ASSERT_TRUE(osh_transport_ion_run_range(&ctx, nn, nn, nn, nn, NULL, 5u, 5u, &completed) == OSH_EINVAL);

    /* nstat == 0. */
    ctx.params.nstat = 0u;
    ASSERT_TRUE(osh_transport_ion_run_range(&ctx, nn, nn, nn, nn, NULL, 0u, 10u, &completed) == OSH_EINVAL);

    /* Out-of-range DELTAE. */
    ctx.params.nstat = 10u;
    ctx.params.deltae = 0.0f;
    ASSERT_TRUE(osh_transport_ion_run_range(&ctx, nn, nn, nn, nn, NULL, 0u, 10u, &completed) == OSH_EINVAL);
}

static double run_boundary_step_energy_loss(int pdg, double energy_mev_per_u) {
    struct osh_particle_pool pool;
    struct particle part;
    struct osh_zone_ref zone_ref;
    struct osh_step_segment segment;
    struct osh_transport_context ctx;
    struct osh_material_runtime mat_rt;
    struct osh_scoring_runtime score_rt;
    struct osh_rng rng;
    float rho_f[3];                /* Runtime density table [g/cm3]. */
    float zero_material[3];        /* Unused material scalar table entries. */
    float range_csda[12];          /* [3 materials][projectile H,He][energy 0,1]. */
    unsigned int projectile_z[2];  /* Runtime projectile columns: H and He. */
    unsigned int projectile_a[2];  /* Representative isotope A for each Z column. */
    double projectile_mass_mev[2]; /* Representative masses; part.mass should override. */
    double e0;                     /* Entry total kinetic energy [MeV]. */

    ASSERT_TRUE(osh_particle_from_pdg(&part, pdg) == 1);
    ASSERT_TRUE(part.a > 0u);
    e0 = energy_mev_per_u * (double) part.a;

    projectile_z[0] = 1u;
    projectile_z[1] = 2u;
    projectile_a[0] = 1u;
    projectile_a[1] = 4u;
    projectile_mass_mev[0] = 938.272046;
    projectile_mass_mev[1] = 3727.379378;

    memset(range_csda, 0, sizeof(range_csda));
    range_csda[8] = 0.0f;
    range_csda[9] = 100.0f;
    range_csda[10] = 0.0f;
    range_csda[11] = 400.0f;

    rho_f[0] = 0.0f;
    rho_f[1] = 0.0f;
    rho_f[2] = 1.19f;
    zero_material[0] = 0.0f;
    zero_material[1] = 0.0f;
    zero_material[2] = 0.0f;

    memset(&mat_rt, 0, sizeof(mat_rt));
    mat_rt.emin = 1.0e-6;
    mat_rt.emax = energy_mev_per_u;
    mat_rt.log_emin = log(mat_rt.emin);
    mat_rt.inv_dlog = 1.0 / log(mat_rt.emax / mat_rt.emin);
    mat_rt.nmaterials = 3u;
    mat_rt.nprojectiles = 2u;
    mat_rt.nenergy = 2u;
    mat_rt.range_csda = range_csda;
    mat_rt.rho = rho_f;
    mat_rt.z_mean = zero_material;
    mat_rt.z_over_a = zero_material;
    mat_rt.rad_length = zero_material;
    mat_rt.moliere_chic2 = zero_material;
    mat_rt.moliere_screen_z = zero_material;
    mat_rt.projectile_z = projectile_z;
    mat_rt.projectile_a = projectile_a;
    mat_rt.projectile_mass_mev = projectile_mass_mev;

    memset(&ctx, 0, sizeof(ctx));
    ctx.params.deltae = 0.99f;
    ctx.params.mcs_mode = OSH_TRANSPORT_MCS_OFF;
    ctx.params.straggling_mode = OSH_TRANSPORT_STRAGGLING_OFF;
    ctx.params.nuclear_inelastic = 0;
    ctx.params.nuclear_elastic = 0;

    memset(&zone_ref, 0, sizeof(zone_ref));
    zone_ref.zone_idx = 0u;
    zone_ref.material_idx = 2u;

    segment.ds = 10.0;

    memset(&score_rt, 0, sizeof(score_rt));
    osh_rng_init(&rng, OSH_RNG_TYPE_PCG32, 123u, 0u);

    ASSERT_TRUE(osh_particle_pool_init(&pool, 1u) == OSH_OK);
    pool.n = 1u;
    pool.x[0] = 0.0;
    pool.y[0] = 0.0;
    pool.z[0] = 0.0;
    pool.ux[0] = 0.0;
    pool.uy[0] = 0.0;
    pool.uz[0] = 1.0;
    pool.e[0] = e0;
    pool.wt[0] = 1.0;
    pool.prim_idx[0] = 0u;
    pool.gen[0] = 0u;
    pool.species[0] = &part;

    ASSERT_TRUE(
        osh_transport_ion_step(&pool, 0u, &zone_ref, &segment, 1u, NULL, &ctx, NULL, &mat_rt, &score_rt, NULL, &rng)
        == OSH_OK);
    ASSERT_TRUE(pool.e[0] > 0.0);

    e0 -= pool.e[0];
    osh_particle_pool_free(&pool);
    return e0;
}

static void test_nondefault_isotope_range_scaling(void) {
    double proton_loss;
    double deuteron_loss;
    double he3_loss;
    double he4_loss;

    proton_loss = run_boundary_step_energy_loss(OSH_PART_PDG_PROTON, 100.0);
    deuteron_loss = run_boundary_step_energy_loss(OSH_PART_PDG_DEUTERON, 100.0);
    he3_loss = run_boundary_step_energy_loss(OSH_PART_PDG_HE3, 100.0);
    he4_loss = run_boundary_step_energy_loss(OSH_PART_PDG_HE4, 100.0);

    ASSERT_TRUE(fabs(deuteron_loss - proton_loss) < 1.0e-4);
    ASSERT_TRUE(fabs(he3_loss - he4_loss) < 1.0e-4);
}

int main(void) {
    test_transport_run_rejects_bad_args();
    test_ion_run_range_validates_arguments();
    test_nondefault_isotope_range_scaling();
    return 0;
}
