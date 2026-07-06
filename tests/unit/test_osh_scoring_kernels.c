#include <math.h>

#include "scoring/runtime/osh_scoring_kernels.h"
#include "test_assert.h"

/* The kernels are exact one-line expressions; compare bit-exactly (fabs<=0 avoids
 * a -Wfloat-equal on ==) against the same expression to lock the op order. */
#define ASSERT_EQ(a, b) ASSERT_TRUE(fabs((a) - (b)) <= 0.0)

static void test_step_kernels(void) {
    /* Energy along a step: fraction path/score_len of de. */
    ASSERT_EQ(osh_kernel_step_energy(2.0, 3.5, 8.0), 2.0 * (3.5 / 8.0));
    /* Fluence: the raw track length (÷volume happens in postprocess). */
    ASSERT_EQ(osh_kernel_step_fluence(3.5), 3.5);
    /* Dose: path length × the pre-gathered dose_scale. */
    ASSERT_EQ(osh_kernel_step_dose(3.5, 0.25), 3.5 * 0.25);
    /* Two-pass weights for LET/Qeff. */
    ASSERT_EQ(osh_kernel_dose_weight(2.0, 3.5, 8.0), 2.0 * 3.5 / 8.0);
    ASSERT_EQ(osh_kernel_track_weight(4.0, 3.5, 8.0), 4.0 * 3.5 / 8.0);
}

static void test_point_kernels(void) {
    /* Point energy: the whole release, no path fraction. */
    ASSERT_EQ(osh_kernel_point_energy(2.0), 2.0);
    /* Point dose: energy per mass, × dose-to-medium ratio (1.0 = no override). */
    ASSERT_EQ(osh_kernel_point_dose(2.0, 1.0, 1.0), 2.0);
    ASSERT_EQ(osh_kernel_point_dose(4.0, 2.0, 1.5), (4.0 / 2.0) * 1.5);
}

/* Step and point are genuinely different functions, not one via a path=1 hack. */
static void test_step_and_point_differ(void) {
    double const de = 2.0;

    /* Along a partial step the energy kernel carries the path fraction... */
    ASSERT_TRUE(osh_kernel_step_energy(de, 3.5, 8.0) < osh_kernel_point_energy(de));
    /* ...and they coincide only in the degenerate point limit (path == score_len),
     * which is exactly why the old path=1 reuse was bit-identical. */
    ASSERT_EQ(osh_kernel_step_energy(de, 1.0, 1.0), osh_kernel_point_energy(de));
    /* Dose: point is de/rho (no path); step needs a path and a gathered scale. */
    ASSERT_EQ(osh_kernel_point_dose(de, 4.0, 1.0), de / 4.0);
}

int main(void) {
    test_step_kernels();
    test_point_kernels();
    test_step_and_point_differ();
    return 0;
}
