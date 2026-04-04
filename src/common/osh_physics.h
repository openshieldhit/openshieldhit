#ifndef OSH_PHYSICS_H
#define OSH_PHYSICS_H

#include <math.h>

/* Relativistic kinematics helpers.
 *
 * All functions operate in natural units:
 *   energy / mass  [MeV or MeV/c²]
 *   momentum       [MeV/c]
 *
 * Mass must be in MeV/c² — use part->mass_MeV (which is A * OSH_AMU + binding
 * correction), not A * OSH_AMU directly.  Converting once at particle
 * initialisation and storing MeV avoids per-call multiplications in the hot
 * path.
 *
 * All functions are static inline so the compiler can inline and optimise them
 * at the call site with zero function-call overhead.  The header has no
 * corresponding .c file.
 *
 * Notation used throughout:
 *   T  — total kinetic energy  [MeV]
 *   p  — total momentum        [MeV/c]
 *   m  — rest mass             [MeV/c²]
 *   E  — total energy = T + m  [MeV]
 *
 * The fundamental relation is: E² = (pc)² + (mc²)²
 * i.e.  (T + m)² = p² + m²  (c = 1). */

/* Momentum from kinetic energy and rest mass.
 *
 * p = sqrt((T + m)^2 - m^2)
 *   = sqrt(T^2 + 2*T*m)
 *
 * The second form avoids catastrophic cancellation for low-energy particles
 * where T << m.
 *
 * @param tkin  total kinetic energy [MeV], must be >= 0
 * @param mass  rest mass [MeV/c²], must be > 0
 * @return total momentum [MeV/c] */
static inline double osh_physics_momentum(double tkin, double mass) {
    return sqrt(tkin * (tkin + 2.0 * mass));
}

/* Kinetic energy from momentum and rest mass.
 *
 * T = sqrt(p^2 + m^2) - m
 *
 * @param mom   total momentum [MeV/c], must be >= 0
 * @param mass  rest mass [MeV/c²], must be > 0
 * @return total kinetic energy [MeV] */
static inline double osh_physics_tkin(double mom, double mass) {
    return sqrt(mom * mom + mass * mass) - mass;
}

/* Total energy (kinetic + rest mass).
 *
 * E = T + m
 *
 * @param tkin  total kinetic energy [MeV]
 * @param mass  rest mass [MeV/c²]
 * @return total energy [MeV] */
static inline double osh_physics_total_energy(double tkin, double mass) {
    return tkin + mass;
}

/* Lorentz factor gamma.
 *
 * gamma = E / m = (T + m) / m = 1 + T/m
 *
 * @param tkin  total kinetic energy [MeV], must be >= 0
 * @param mass  rest mass [MeV/c²], must be > 0
 * @return dimensionless Lorentz factor (>= 1) */
static inline double osh_physics_gamma(double tkin, double mass) {
    return 1.0 + tkin / mass;
}

/* Velocity as fraction of speed of light, beta = v/c.
 *
 * beta = p / E = p / (T + m)
 *
 * @param tkin  total kinetic energy [MeV], must be >= 0
 * @param mass  rest mass [MeV/c²], must be > 0
 * @return beta in [0, 1) */
static inline double osh_physics_beta(double tkin, double mass) {
    double p = osh_physics_momentum(tkin, mass);
    return p / (tkin + mass);
}

#endif /* OSH_PHYSICS_H */
