#include "common/osh_sysinfo.h"

#include <stdio.h>

/*
 * Platform resource detection.  Each branch is self-contained; the goal is a
 * best-effort snapshot, never a hard failure.  Unknown values stay 0.
 *
 * Memory semantics we aim for:
 *   ram_total_bytes     = usable RAM ceiling: physical RAM, but capped by a
 *                         cgroup/container limit (Linux) or the heap maximum
 *                         (WASM) when those are smaller — because that is the
 *                         real ceiling a run will hit.
 *   ram_available_bytes = what could be allocated right now (free plus easily
 *                         reclaimable memory), as reported by the OS.
 */

#if defined(__EMSCRIPTEN__)
/* ---- WebAssembly (Emscripten) ------------------------------------------- */
#include <emscripten/emscripten.h>
#include <emscripten/threading.h>

static void sysinfo_query_platform(struct osh_sysinfo *out) {
    /* The wasm linear memory maximum is the hard ceiling (wasm32 also caps the
     * address space at ~4 GiB regardless of host RAM). */
    size_t const heap_max = emscripten_get_heap_max();
    size_t const heap_used = emscripten_get_heap_size();

    int const cores = emscripten_num_logical_cores();
    out->logical_cores = (cores > 0) ? (unsigned) cores : 0u;
    out->ram_total_bytes = (uint64_t) heap_max;
    out->ram_available_bytes = (heap_max > heap_used) ? (uint64_t) (heap_max - heap_used) : 0u;
}

#elif defined(_WIN32)
/* ---- Windows ------------------------------------------------------------- */
#include <windows.h>

static void sysinfo_query_platform(struct osh_sysinfo *out) {
    MEMORYSTATUSEX mem;
    DWORD cores;

    cores = GetActiveProcessorCount(ALL_PROCESSOR_GROUPS);
    out->logical_cores = (cores > 0) ? (unsigned) cores : 0u;

    mem.dwLength = sizeof(mem);
    if (GlobalMemoryStatusEx(&mem)) {
        out->ram_total_bytes = (uint64_t) mem.ullTotalPhys;
        out->ram_available_bytes = (uint64_t) mem.ullAvailPhys;
    }
}

#elif defined(__APPLE__)
/* ---- macOS --------------------------------------------------------------- */
#include <mach/mach.h>
#include <mach/mach_host.h>
#include <sys/sysctl.h>
#include <unistd.h>

static void sysinfo_query_platform(struct osh_sysinfo *out) {
    int ncpu = 0;
    size_t len = sizeof(ncpu);
    uint64_t memsize = 0;

    if (sysctlbyname("hw.logicalcpu", &ncpu, &len, NULL, 0) == 0 && ncpu > 0) {
        out->logical_cores = (unsigned) ncpu;
    }

    len = sizeof(memsize);
    if (sysctlbyname("hw.memsize", &memsize, &len, NULL, 0) == 0) {
        out->ram_total_bytes = memsize;
    }

    /* Available ≈ (free + inactive) pages: inactive pages are reclaimable. */
    {
        mach_port_t host = mach_host_self();
        vm_size_t page_size = 0;
        vm_statistics64_data_t vm;
        mach_msg_type_number_t count = HOST_VM_INFO64_COUNT;

        if (host_page_size(host, &page_size) == KERN_SUCCESS
            && host_statistics64(host, HOST_VM_INFO64, (host_info64_t) &vm, &count) == KERN_SUCCESS) {
            uint64_t const reclaimable = (uint64_t) vm.free_count + (uint64_t) vm.inactive_count;
            out->ram_available_bytes = reclaimable * (uint64_t) page_size;
        }
    }
}

#elif defined(__linux__)
/* ---- Linux --------------------------------------------------------------- */
#include <unistd.h>

/* Read a single unsigned integer from a one-line sysfs/proc file.
 * Returns 1 and sets *out on success; 0 if the file is absent or holds a
 * non-numeric sentinel such as cgroup-v2's "max". */
static int read_u64_file(char const *path, uint64_t *out) {
    FILE *f = fopen(path, "r");
    unsigned long long v;
    int got;

    if (!f) {
        return 0;
    }
    got = fscanf(f, "%llu", &v);
    fclose(f);
    if (got != 1) {
        return 0; /* e.g. "max" → unlimited */
    }
    *out = (uint64_t) v;
    return 1;
}

/* Parse MemAvailable (kB) from /proc/meminfo; returns bytes or 0 if absent. */
static uint64_t linux_mem_available(void) {
    FILE *f = fopen("/proc/meminfo", "r");
    char line[256];
    uint64_t result = 0u;

    if (!f) {
        return 0u;
    }
    while (fgets(line, sizeof(line), f)) {
        unsigned long long kb;
        if (sscanf(line, "MemAvailable: %llu kB", &kb) == 1) {
            result = (uint64_t) kb * 1024u;
            break;
        }
    }
    fclose(f);
    return result;
}

/* A cgroup limit this large is the kernel's "unlimited" sentinel, not a real
 * cap; ignore anything at/above a generous threshold. */
#define OSH_CGROUP_UNLIMITED_MIN ((uint64_t) 1 << 62)

static uint64_t linux_cgroup_limit(void) {
    uint64_t v = 0u;

    /* cgroup v2 first, then v1. */
    if (read_u64_file("/sys/fs/cgroup/memory.max", &v) && v < OSH_CGROUP_UNLIMITED_MIN) {
        return v;
    }
    if (read_u64_file("/sys/fs/cgroup/memory/memory.limit_in_bytes", &v) && v < OSH_CGROUP_UNLIMITED_MIN) {
        return v;
    }
    return 0u; /* unlimited / not present */
}

static void sysinfo_query_platform(struct osh_sysinfo *out) {
    long const ncpu = sysconf(_SC_NPROCESSORS_ONLN);
    long const pages = sysconf(_SC_PHYS_PAGES);
    long const page_size = sysconf(_SC_PAGESIZE);
    uint64_t physical = 0u;
    uint64_t cgroup;
    uint64_t avail;

    if (ncpu > 0) {
        out->logical_cores = (unsigned) ncpu;
    }
    if (pages > 0 && page_size > 0) {
        physical = (uint64_t) pages * (uint64_t) page_size;
    }

    /* Honor a container/cgroup limit when it is smaller than physical RAM —
     * the JVM's well-known container-OOM bug was exactly ignoring this. */
    cgroup = linux_cgroup_limit();
    if (physical == 0u) {
        out->ram_total_bytes = cgroup;
    } else if (cgroup > 0u && cgroup < physical) {
        out->ram_total_bytes = cgroup;
    } else {
        out->ram_total_bytes = physical;
    }

    avail = linux_mem_available();
    /* Never report "available" above the (possibly cgroup-capped) total. */
    if (avail > 0u && out->ram_total_bytes > 0u && avail > out->ram_total_bytes) {
        avail = out->ram_total_bytes;
    }
    out->ram_available_bytes = avail;
}

#else
/* ---- Unknown platform: best-effort zero fill ----------------------------- */
static void sysinfo_query_platform(struct osh_sysinfo *out) {
    (void) out; /* all fields already zeroed by the caller */
}
#endif

enum osh_status osh_sysinfo_query(struct osh_sysinfo *out) {
    if (!out) {
        return OSH_EINVAL;
    }

    out->logical_cores = 0u;
    out->ram_total_bytes = 0u;
    out->ram_available_bytes = 0u;
    out->gpu_count = 0u;

    sysinfo_query_platform(out);
    return OSH_OK;
}

void osh_sysinfo_format_bytes(uint64_t bytes, char *buf, size_t buflen) {
    static char const *const units[] = {"B", "KiB", "MiB", "GiB", "TiB", "PiB"};
    size_t const n_units = sizeof(units) / sizeof(units[0]);
    double value = (double) bytes;
    size_t unit = 0u;

    if (!buf || buflen == 0u) {
        return;
    }

    while (value >= 1024.0 && unit + 1u < n_units) {
        value /= 1024.0;
        ++unit;
    }

    if (unit == 0u) {
        /* Exact byte count — no decimals. */
        (void) snprintf(buf, buflen, "%llu %s", (unsigned long long) bytes, units[unit]);
    } else {
        (void) snprintf(buf, buflen, "%.1f %s", value, units[unit]);
    }
}
