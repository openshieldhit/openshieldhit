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
 * @enum osh_rng_purpose
 *
 * @brief Independent sub-stream selector for per-history seeding.
 *
 * @details
 * A single history index can drive several mutually independent streams by
 * mixing a distinct @ref osh_rng_purpose into the stream id (see
 * osh_rng_seed_history()).  This keeps source sampling reproducible and
 * decoupled from stochastic transport: the same primary phase space is drawn
 * regardless of which physics options (MSCAT/STRAGG/NUCRE) are enabled.
 */
enum osh_rng_purpose {
    OSH_RNG_PURPOSE_BEAM = 0,    /**< Source / beam phase-space sampling. */
    OSH_RNG_PURPOSE_PHYSICS = 1, /**< In-transport stochastic physics. */
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
 * @brief Run-wide RNG seeding context for per-history streams.
 *
 * @details
 * Carries the two values needed to derive an independent stream for any
 * history index: the engine @p type and the run @p seed.  @p seed already
 * combines RNDSEED and RNDOFFSET (see @ref osh_rng_seeding_init), so distinct
 * RNDOFFSET values from the same RNDSEED select decorrelated stream families
 * over the same history-index range.  Process/MPI/worker splitting is a
 * separate, orthogonal concern: callers pass disjoint history-index ranges
 * explicitly (e.g. @ref osh_beam_runtime_fill_pool_at's global_prim_base), so
 * rank r owns history indices [r·N, (r+1)·N) regardless of @p seed.
 */
struct osh_rng_seeding {
    enum osh_rng_type type; /**< Engine used for every stream in the run. */
    uint64_t seed;          /**< Run seed; already combines RNDSEED and RNDOFFSET (@ref osh_rng_seeding_init). */
};

/**
 * @brief Build a seeding context from RNDSEED and RNDOFFSET.
 *
 * @details
 * @p rndoffset == 0 (no `-N`/`--seedoffset`, or `-N 0`) is the identity case:
 * @p seed comes out equal to @p rndseed, so the ordinary single-run path is
 * unchanged. A non-zero @p rndoffset instead hashes the pair through the same
 * SplitMix64-style mixer used for per-history streams (`rng_mix_stream()` in
 * osh_rng.c), rather than adding it to @p rndseed. Addition would let two
 * independently-chosen configurations collide: (RNDSEED=S, RNDOFFSET=k) and
 * (RNDSEED=S+k, RNDOFFSET=0) sum to the same value and would therefore
 * produce byte-identical output, silently merging as though the two runs
 * were independent replicas and double-counting histories. The hash-mix
 * reduces that to a generic 64-bit hash coincidence (~2^-64), the same
 * residual risk already accepted for every other stream-separation axis in
 * this module.
 *
 * @param[out] seeding   Seeding context to initialise.
 * @param[in]  type      Engine to use for every stream in the run.
 * @param[in]  rndseed   Base RNG seed (RNDSEED).
 * @param[in]  rndoffset Independent-stream selector (RNDOFFSET / `-N`).
 */
void osh_rng_seeding_init(struct osh_rng_seeding *seeding,
                          enum osh_rng_type type,
                          uint64_t rndseed,
                          uint64_t rndoffset);

/**
 * @brief Seed an RNG for one history, keyed by its global index and purpose.
 *
 * @details
 * Derives a stream id by mixing (@p seed, @p hist_index, @p purpose) through a
 * SplitMix64 finaliser, then initialises @p rng on that stream.  The resulting
 * stream depends only on its key, never on execution order, so the same
 * history sees the same draws regardless of pool capacity, thread, or rank.
 * Distinct @ref osh_rng_purpose values yield independent streams for the same
 * history.
 *
 * @param rng        RNG state to initialise.
 * @param type       Engine type.
 * @param seed       Run seed.
 * @param hist_index Global history index (e.g. hist_lo + worker-local index).
 * @param purpose    Sub-stream selector (see @ref osh_rng_purpose).
 */
void osh_rng_seed_history(
    struct osh_rng *rng, enum osh_rng_type type, uint64_t seed, uint64_t hist_index, enum osh_rng_purpose purpose);

/**
 * @brief Derive an independent child stream from a parent, keyed by ordinal.
 *
 * @details
 * Seeds @p child by hashing @p parent's current internal state, keyed by the
 * child's @p ordinal.  Crucially it does **not** consume a draw from @p parent:
 * splitting reads the parent's state but never advances it.  A secondary that
 * is later dropped, reordered, or lost to pool overflow therefore cannot shift
 * its parent's or its siblings' streams, so reproducibility is independent of
 * pool occupancy and wavefront scheduling (issue #213; design in #148).
 *
 * The child stream is a pure function of the parent's current state (itself a
 * pure function of the parent's lineage) and @p ordinal.  The lineage key is
 * hashed from the parent's raw state words, which the engine permutes before
 * emitting, so a child seed is no longer drawn from — or structurally
 * correlated with — the parent's own subsequent output (issue #299).  Siblings
 * produced by one event are separated by giving each a distinct @p ordinal (its
 * index in the event's secondary list).
 *
 * @param child   RNG state to initialise (engine type inherited from parent).
 * @param parent  Parent RNG; read only, never advanced.
 * @param ordinal Zero-based index of this child among its siblings.
 */
void osh_rng_split(struct osh_rng *child, struct osh_rng const *parent, uint64_t ordinal);

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
 * @brief Sample a Poisson-distributed integer with mean lambda (Knuth algorithm).
 *
 * @details
 * Uses Knuth's product-of-uniforms method, O(lambda) on average.
 * Suitable for modest lambda values. The implementation has a finite loop
 * guard for pathological inputs; callers that need a hard output cap should
 * clamp the returned value themselves.
 *
 * @param rng     RNG state (mutated).
 * @param lambda  Mean of the Poisson distribution; must be >= 0.
 * @returns       Non-negative integer drawn from Poisson(lambda).
 */
int osh_rng_poisson(struct osh_rng *rng, double lambda);

/**
 * @brief Convenience alias for generating a double in the range [0, 1).
 *
 * @param rng Pointer to the RNG state.
 *
 * @return Double in the range [0, 1).
 */
static inline double osh_rng(struct osh_rng *rng) {
    return osh_rng_double(rng);
}

#endif /* OSH_RNG_H */
