#ifndef OSH_RNG_H
#define OSH_RNG_H

/**
 * @file osh_rng.h
 * @brief OpenShieldHIT RNG (engine + distributions)
 *
 * Design goals:
 * - Stack-only state (no heap allocation, no pointers required)
 * - Runtime engine selection (switch-based dispatch in implementation)
 * - Fast uniform draws (u32/u64/f32/f64)
 * - Fast Gaussian sampling (Box-Muller with cached spare)
 *
 * Notes:
 * - "seed" selects the run; "stream" (a.k.a. sequence id) selects an
 *   independent random sequence for parallelism (thread/history lanes).
 */

#include <stdint.h>

/** Forward declaration for engine API prototypes below. */
struct osh_rng;

/**
 * @name Engine-specific functions
 * @{
 */

/**
 * @brief Initialize PCG32 engine.
 *
 * @param rng Pointer to the RNG state.
 * @param seed Seed value.
 * @param stream Stream/sequence ID.
 */
void osh_rng_pcg32_init(struct osh_rng *rng, uint64_t seed, uint64_t stream);

/**
 * @brief Generate a 32-bit unsigned integer using PCG32 engine.
 *
 * @param rng Pointer to the RNG state.
 *
 * @return 32-bit unsigned integer.
 */
uint32_t osh_rng_pcg32_u32(struct osh_rng *rng);

/**
 * @brief Initialize xoshiro256** engine.
 *
 * @param rng Pointer to the RNG state.
 * @param seed Seed value.
 * @param stream Stream/sequence ID.
 */
void osh_rng_xoshiro256ss_init(struct osh_rng *rng, uint64_t seed, uint64_t stream);

/**
 * @brief Generate a 64-bit unsigned integer using xoshiro256** engine.
 *
 * @param rng Pointer to the RNG state.
 *
 * @return 64-bit unsigned integer.
 */
uint64_t osh_rng_xoshiro256ss_u64(struct osh_rng *rng);

/** @} */

/**
 * @enum osh_rng_type
 *
 * @brief Enumeration of RNG engine types.
 */
enum osh_rng_type {
    OSH_RNG_TYPE_PCG32 = 1,        /**< PCG32 engine */
    OSH_RNG_TYPE_XOSHIRO256SS = 2, /**< xoshiro256** engine */
};

/**
 * @struct osh_rng
 *
 * @brief RNG state container.
 *
 * Keep this on the stack or embed in other state objects.
 *
 * The "gauss_has_spare" / "gauss_spare" cache is used by osh_rng_gauss*()
 * to return two normal variates per underlying transform.
 */
struct osh_rng {
    enum osh_rng_type type; /**< RNG engine type */

    union {
        struct {
            uint64_t state; /**< RNG state */
            uint64_t inc;   /**< Increment value (must be odd) */
        } pcg32;

        struct {
            uint64_t s[4]; /**< RNG state array */
        } xoshiro256ss;
    } u;

    double gauss_spare;  /**< Cached spare value for Gaussian sampling */
    int gauss_has_spare; /**< Flag indicating if spare value is available */
};

/**
 * @brief Initialize RNG with selected engine, seed, and stream/sequence ID.
 *
 * @param rng Pointer to the RNG state.
 * @param type RNG engine type.
 * @param seed Seed value.
 * @param stream Stream/sequence ID.
 */
void osh_rng_init(struct osh_rng *rng, enum osh_rng_type type, uint64_t seed, uint64_t stream);

/**
 * @brief Generate a 32-bit unsigned integer.
 *
 * @param rng Pointer to the RNG state.
 *
 * @return 32-bit unsigned integer.
 */
uint32_t osh_rng_u32(struct osh_rng *rng);

/**
 * @brief Generate a 64-bit unsigned integer.
 *
 * @param rng Pointer to the RNG state.
 *
 * @return 64-bit unsigned integer.
 */
uint64_t osh_rng_u64(struct osh_rng *rng);

/**
 * @brief Generate a float in the range [0, 1).
 *
 * @param rng Pointer to the RNG state.
 *
 * @return Float in the range [0, 1).
 */
float osh_rng_float(struct osh_rng *rng);

/**
 * @brief Generate a double in the range [0, 1).
 *
 * @param rng Pointer to the RNG state.
 *
 * @return Double in the range [0, 1).
 */
double osh_rng_double(struct osh_rng *rng);

/**
 * @brief Generate a standard normal random variable (N(0,1)).
 *
 * @param rng Pointer to the RNG state.
 *
 * @return Standard normal random variable.
 */
double osh_rng_gauss01(struct osh_rng *rng);

/**
 * @brief Generate a normal random variable (N(mu, sigma)).
 *
 * @param rng Pointer to the RNG state.
 * @param mu Mean of the distribution.
 * @param sigma Standard deviation of the distribution.
 *
 * @return Normal random variable.
 */
double osh_rng_gauss(struct osh_rng *rng, double mu, double sigma);

/**
 * @brief Generate an array of doubles in the range [0, 1).
 *
 * @param rng Pointer to the RNG state.
 * @param x Pointer to the output array.
 * @param n Number of elements to generate.
 */
void osh_rng_double_vec(struct osh_rng *rng, double *restrict x, int n);

/**
 * @brief Generate an array of floats in the range [0, 1).
 *
 * @param rng Pointer to the RNG state.
 * @param x Pointer to the output array.
 * @param n Number of elements to generate.
 */
void osh_rng_float_vec(struct osh_rng *rng, float *restrict x, int n);

/**
 * @brief Generate an array of 32-bit unsigned integers.
 *
 * @param rng Pointer to the RNG state.
 * @param x Pointer to the output array.
 * @param n Number of elements to generate.
 */
void osh_rng_u32_vec(struct osh_rng *rng, uint32_t *restrict x, int n);

/**
 * @brief Generate an array of standard normal random variables (N(0,1)).
 *
 * @param rng Pointer to the RNG state.
 * @param x Pointer to the output array.
 * @param n Number of elements to generate.
 */
void osh_rng_gauss01_vec(struct osh_rng *rng, double *restrict x, int n);

/**
 * @brief Generate an array of normal random variables (N(mu, sigma)).
 *
 * @param rng Pointer to the RNG state.
 * @param mu Mean of the distribution.
 * @param sigma Standard deviation of the distribution.
 * @param x Pointer to the output array.
 * @param n Number of elements to generate.
 */
void osh_rng_gauss_vec(struct osh_rng *rng, double mu, double sigma, double *restrict x, int n);

/**
 * @brief Convenience alias for generating a double in the range [0, 1).
 *
 * @param rng Pointer to the RNG state.
 *
 * @return Double in the range [0, 1).
 */
static inline double osh_rng(struct osh_rng *rng) { return osh_rng_double(rng); }

#endif /* OSH_RNG_H */
