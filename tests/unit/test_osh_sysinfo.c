/*
 * Unit tests for osh_sysinfo: best-effort host resource detection and the
 * byte-formatting helper.  Detection values are platform-dependent, so the
 * assertions check invariants (cores >= 1, total > 0 on supported desktop
 * OSes, available <= total) rather than exact numbers.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_sysinfo.h"
#include "openshieldhit/status.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

#define ASSERT_STREQ(a, b)                                                                                             \
    do {                                                                                                               \
        char const *_a = (a);                                                                                          \
        char const *_b = (b);                                                                                          \
        if (strcmp(_a, _b) != 0) {                                                                                     \
            fprintf(stderr, "ASSERT FAILED: \"%s\" != \"%s\" (%s:%d)\n", _a, _b, __FILE__, __LINE__);                  \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

static void test_query_null_is_einval(void) {
    ASSERT_TRUE(osh_sysinfo_query(NULL) == OSH_EINVAL);
}

static void test_query_invariants(void) {
    struct osh_sysinfo info;

    /* Sentinel so we can detect that the query zero-initialises before filling. */
    info.gpu_count = 12345u;

    ASSERT_TRUE(osh_sysinfo_query(&info) == OSH_OK);

    /* GPU detection is not implemented yet: must be reset to 0. */
    ASSERT_TRUE(info.gpu_count == 0u);

    /* On every CI target (Linux/macOS/Windows) these are reported. */
    ASSERT_TRUE(info.logical_cores >= 1u);
    ASSERT_TRUE(info.ram_total_bytes > 0u);

    /* Available is best-effort; when known it cannot exceed total. */
    if (info.ram_available_bytes > 0u) {
        ASSERT_TRUE(info.ram_available_bytes <= info.ram_total_bytes);
    }
}

static void test_format_bytes(void) {
    char buf[32];

    osh_sysinfo_format_bytes(0u, buf, sizeof(buf));
    ASSERT_STREQ(buf, "0 B");

    osh_sysinfo_format_bytes(512u, buf, sizeof(buf));
    ASSERT_STREQ(buf, "512 B");

    osh_sysinfo_format_bytes(1023u, buf, sizeof(buf));
    ASSERT_STREQ(buf, "1023 B");

    osh_sysinfo_format_bytes(1024u, buf, sizeof(buf));
    ASSERT_STREQ(buf, "1.0 KiB");

    osh_sysinfo_format_bytes(1536u, buf, sizeof(buf));
    ASSERT_STREQ(buf, "1.5 KiB");

    osh_sysinfo_format_bytes((uint64_t) 1u << 20, buf, sizeof(buf));
    ASSERT_STREQ(buf, "1.0 MiB");

    osh_sysinfo_format_bytes((uint64_t) 3u << 30, buf, sizeof(buf));
    ASSERT_STREQ(buf, "3.0 GiB");

    osh_sysinfo_format_bytes((uint64_t) 1u << 40, buf, sizeof(buf));
    ASSERT_STREQ(buf, "1.0 TiB");
}

static void test_format_bytes_no_overflow_on_tiny_buffer(void) {
    char buf[4];

    /* Must NUL-terminate and not write past the buffer; content is truncated. */
    osh_sysinfo_format_bytes((uint64_t) 3u << 30, buf, sizeof(buf));
    ASSERT_TRUE(buf[sizeof(buf) - 1] == '\0');

    /* buflen == 0 must be a no-op (no write, no crash). */
    osh_sysinfo_format_bytes(123u, buf, 0u);
}

int main(void) {
    test_query_null_is_einval();
    test_query_invariants();
    test_format_bytes();
    test_format_bytes_no_overflow_on_tiny_buffer();
    return 0;
}
