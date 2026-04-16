#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "openshieldhit/version.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

static void test_version_string_defined(void) {
    char const *v = osh_version_string();
    ASSERT_TRUE(v != NULL);
    ASSERT_TRUE(strlen(v) > 0);
}

static void test_version_components_non_negative(void) {
    ASSERT_TRUE(osh_version_major() >= 0);
    ASSERT_TRUE(osh_version_minor() >= 0);
    ASSERT_TRUE(osh_version_patch() >= 0);
}

static void test_version_components_in_string(void) {
    char expected[64];
    snprintf(expected, sizeof(expected), "%d.%d.%d", osh_version_major(), osh_version_minor(), osh_version_patch());
    ASSERT_TRUE(strstr(osh_version_string(), expected) != NULL);
}

static int run_named_test(char const *name) {
    if (strcmp(name, "version_string_defined") == 0) {
        test_version_string_defined();
        return 0;
    }
    if (strcmp(name, "version_components_non_negative") == 0) {
        test_version_components_non_negative();
        return 0;
    }
    if (strcmp(name, "version_components_in_string") == 0) {
        test_version_components_in_string();
        return 0;
    }
    return 1;
}

int main(int argc, char *argv[]) {
    if (argc == 2) {
        return run_named_test(argv[1]);
    }

    test_version_string_defined();
    test_version_components_non_negative();
    test_version_components_in_string();
    return 0;
}
