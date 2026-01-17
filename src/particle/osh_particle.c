
#include "osh_particle.h"
#include <stdio.h>

/**
 * @brief Print particle struct.
 *
 * @param[in] part - particle to be printed
 *
 * @author Niels Bassler
 */
void osh_particle_print(struct particle part) {
    printf("Print particle:\n");
    printf("  PDG    = %i\n", part.pdg);
    printf("  id     = %i\n", part.id);
    printf("  Z      = %i\n", part.z);
    printf("  A      = %i\n", part.a);
    printf("  amass  = %12.5f MeV/c**2\n", part.amass);
    printf("  amu    = %12.5f u\n", part.amu);
    printf("  gen    = %i\n", part.gen);
    printf("  nprim  = %i\n", part.nprim);
    printf("  weight = %f\n", part.weight);
}
