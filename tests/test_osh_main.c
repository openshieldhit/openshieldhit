#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_version.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

static void test_version_string_defined(void) {
    char const *v = OSH_VERSION;
    ASSERT_TRUE(v != NULL);
    ASSERT_TRUE(strlen(v) > 0);
}

static void test_version_components_non_negative(void) {
    ASSERT_TRUE(OSH_VERSION_MAJOR >= 0);
    ASSERT_TRUE(OSH_VERSION_MINOR >= 0);
    ASSERT_TRUE(OSH_VERSION_PATCH >= 0);
}

static void test_version_components_in_string(void) {
    char expected[64];
    snprintf(expected, sizeof(expected), "%d.%d.%d", OSH_VERSION_MAJOR, OSH_VERSION_MINOR, OSH_VERSION_PATCH);
    ASSERT_TRUE(strstr(OSH_VERSION, expected) != NULL);
}

int main(void) {
    test_version_string_defined();
    test_version_components_non_negative();
    test_version_components_in_string();
    return 0;
}
