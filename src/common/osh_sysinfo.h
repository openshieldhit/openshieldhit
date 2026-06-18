#ifndef OSH_SYSINFO_H
#define OSH_SYSINFO_H

/**
 * @file osh_sysinfo.h
 * @brief Best-effort host resource detection (CPU cores, RAM; later: GPU).
 *
 * @details
 * This is a thin platform-abstraction utility, the resource-detection
 * counterpart of osh_time.c.  It is the ONLY place in the project that queries
 * the operating system for hardware/resource information, so all OS-specific
 * code (`#ifdef _WIN32 / __APPLE__ / __linux__ / __EMSCRIPTEN__`) is confined
 * here and the rest of the tree stays portable.
 *
 * Design intent (parallelization preparation, see issue #152):
 *   - The CORE library never calls this.  Resource detection is a policy
 *     concern owned by the application (or a library embedder), which may also
 *     inject its own numbers instead of probing the OS.  The core simulation
 *     only ever consumes a plain memory-budget number.
 *   - The query is best-effort: it fills whatever the platform can report and
 *     leaves unknown fields as 0.  It never fails hard.
 *
 * This module currently powers the `--print-resources` CLI action.  A later
 * change (component B of #152) will use the same struct to size a memory
 * budget and gate runs that would otherwise risk out-of-memory.
 */

#include <stddef.h>
#include <stdint.h>

#include "openshieldhit/status.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Snapshot of host resources relevant to sizing a run.
 *
 * @details
 * Any field that the platform cannot report is left as 0 ("unknown"), so
 * callers must treat 0 as "no information" rather than "zero cores / zero
 * bytes".  Memory figures are physical RAM, except that on constrained
 * environments they reflect the binding limit: a Linux cgroup/container limit
 * (taken as the min with physical RAM) or, under WebAssembly, the configured
 * heap maximum.
 */
struct osh_sysinfo {
    unsigned logical_cores;       /**< Online logical CPUs; 0 = unknown. */
    uint64_t ram_total_bytes;     /**< Total usable RAM (physical or cgroup/WASM cap); 0 = unknown. */
    uint64_t ram_available_bytes; /**< Currently available RAM (free + reclaimable); 0 = unknown. */

    /* Reserved for future GPU detection (component of the GPU work in #148);
     * always 0 today so callers can already branch on it. */
    unsigned gpu_count; /**< Number of detected compute GPUs; 0 = none/unknown. */
};

/**
 * @brief Query host resources into @p out (best-effort).
 *
 * @details
 * Fills @p out with whatever the current platform can report, leaving unknown
 * fields 0.  Implemented per platform (Linux, macOS, Windows, Emscripten) with
 * a zero-filled fallback for any other target.
 *
 * @param[out] out  Receives the snapshot; must not be NULL.
 *
 * @returns OSH_OK on success (including the partial/zero-filled case),
 *          OSH_EINVAL if @p out is NULL.
 */
enum osh_status osh_sysinfo_query(struct osh_sysinfo *out);

/**
 * @brief Format a byte count as a short human-readable string (binary units).
 *
 * @details
 * Produces values such as "0 B", "512 KiB", "61.2 GiB" using 1024-based units
 * (B, KiB, MiB, GiB, TiB, PiB).  Values >= 1 KiB are printed with one decimal;
 * bytes are printed as an integer.  Always NUL-terminates if @p buflen > 0.
 *
 * This lives here (rather than in the reporting code) because the same
 * formatting is reused by the resource report and, later, by the memory-budget
 * report (component B of #152).
 *
 * @param[in]  bytes   Value to format.
 * @param[out] buf     Output buffer.
 * @param[in]  buflen  Size of @p buf in bytes.
 */
void osh_sysinfo_format_bytes(uint64_t bytes, char *buf, size_t buflen);

#ifdef __cplusplus
}
#endif

#endif /* OSH_SYSINFO_H */
