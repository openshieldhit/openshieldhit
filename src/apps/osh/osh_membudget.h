#ifndef OSH_MEMBUDGET_H
#define OSH_MEMBUDGET_H

/**
 * @file osh_membudget.h
 * @brief Memory-budget parsing and the default budgeting policy.
 *
 * @details
 * App-layer helpers that turn detected host resources (@ref osh_sysinfo) and an
 * optional user override string into a concrete byte budget.  Kept separate
 * from osh_run.c so the parsing and policy can be unit-tested in isolation.
 *
 * These are pure functions: they never touch the OS or allocate.  Resource
 * detection lives in osh_sysinfo; consuming the budget (estimating scoring
 * memory and refusing over-budget runs) lives in osh_run.c.
 */

#include <stdint.h>

#include "common/osh_sysinfo.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Parse a memory-budget string into an absolute byte count.
 *
 * @details
 * Accepts:
 *   - a plain byte count: `"1048576"`;
 *   - a size with a binary unit suffix (case-insensitive, 1024-based;
 *     `KB`/`KiB` are treated identically): `"512MB"`, `"4GiB"`, `"2g"`;
 *   - a percentage of @p base: `"80%"`.
 *
 * The value may be fractional for unit/percent forms (`"1.5GB"`, `"62.5%"`).
 * Leading/trailing spaces are tolerated.  Negative, empty, non-numeric, or
 * out-of-range inputs are rejected.
 *
 * @param[in]  text  The budget string (must not be NULL).
 * @param[in]  base  Reference size for the `%` form (typically available RAM).
 * @param[out] out_bytes  Receives the parsed byte count on success.
 *
 * @returns 1 on success, 0 if @p text is malformed or out of range.
 */
int osh_membudget_parse(char const *text, uint64_t base, uint64_t *out_bytes);

/**
 * @brief Compute the default memory budget from detected resources.
 *
 * @details
 * Policy: 80% of currently *available* RAM (falling back to total RAM if
 * available is unknown), but never closer than a reserve floor of
 * `max(256 MiB, 5% of total)` to the limit, so the OS and OpenShieldHIT's
 * other (small) allocations keep breathing room.  Basing on *available* rather
 * than total already accounts for memory other processes hold.
 *
 * @param[in] info  Detected host resources.
 *
 * @returns The budget in bytes, or 0 when no memory figure is known (in which
 *          case the caller should skip budget enforcement).
 */
uint64_t osh_membudget_default(struct osh_sysinfo const *info);

#ifdef __cplusplus
}
#endif

#endif /* OSH_MEMBUDGET_H */
