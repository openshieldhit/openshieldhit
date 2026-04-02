#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "openshieldhit/openshieldhit.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

static void test_version_string_defined(void) {
    char const *v = openshieldhit_version_string();
    ASSERT_TRUE(v != NULL);
    ASSERT_TRUE(strlen(v) > 0);
}

static void test_version_components_non_negative(void) {
    ASSERT_TRUE(openshieldhit_version_major() >= 0);
    ASSERT_TRUE(openshieldhit_version_minor() >= 0);
    ASSERT_TRUE(openshieldhit_version_patch() >= 0);
}

static void test_version_components_in_string(void) {
    char expected[64];
    snprintf(expected,
             sizeof(expected),
             "%d.%d.%d",
             openshieldhit_version_major(),
             openshieldhit_version_minor(),
             openshieldhit_version_patch());
    ASSERT_TRUE(strstr(openshieldhit_version_string(), expected) != NULL);
}

static void test_context_create_destroy(void) {
    openshieldhit_context_t *ctx = openshieldhit_context_create();
    ASSERT_TRUE(ctx != NULL);
    ASSERT_TRUE(openshieldhit_last_error(ctx) != NULL);
    openshieldhit_context_destroy(ctx);
}

static void test_setters_return_ok(void) {
    openshieldhit_context_t *ctx = openshieldhit_context_create();
    ASSERT_TRUE(ctx != NULL);

    ASSERT_TRUE(openshieldhit_context_set_workdir(ctx,     "/tmp")         == OPENSHIELDHIT_STATUS_OK);
    ASSERT_TRUE(openshieldhit_context_set_out_dir(ctx,     "/tmp/out")     == OPENSHIELDHIT_STATUS_OK);
    ASSERT_TRUE(openshieldhit_context_set_geo_path(ctx,    "geo.dat")      == OPENSHIELDHIT_STATUS_OK);
    ASSERT_TRUE(openshieldhit_context_set_beam_path(ctx,   "beam.dat")     == OPENSHIELDHIT_STATUS_OK);
    ASSERT_TRUE(openshieldhit_context_set_mat_path(ctx,    "mat.dat")      == OPENSHIELDHIT_STATUS_OK);
    ASSERT_TRUE(openshieldhit_context_set_detect_path(ctx, "detect.dat")   == OPENSHIELDHIT_STATUS_OK);
    ASSERT_TRUE(openshieldhit_context_set_nstat(ctx,       1000ULL)        == OPENSHIELDHIT_STATUS_OK);
    ASSERT_TRUE(openshieldhit_context_set_log_level(ctx,   1)              == OPENSHIELDHIT_STATUS_OK);
    ASSERT_TRUE(openshieldhit_context_set_run_mode(ctx, OPENSHIELDHIT_RUN_VALIDATE) == OPENSHIELDHIT_STATUS_OK);

    openshieldhit_context_destroy(ctx);
}

static void test_setters_reject_null_context(void) {
    ASSERT_TRUE(openshieldhit_context_set_workdir(NULL, "/tmp")  == OPENSHIELDHIT_STATUS_INVALID_ARGUMENT);
    ASSERT_TRUE(openshieldhit_context_set_nstat(NULL, 100ULL)    == OPENSHIELDHIT_STATUS_INVALID_ARGUMENT);
    ASSERT_TRUE(openshieldhit_context_set_run_mode(NULL, OPENSHIELDHIT_RUN_VALIDATE) == OPENSHIELDHIT_STATUS_INVALID_ARGUMENT);
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
    if (strcmp(name, "context_create_destroy") == 0) {
        test_context_create_destroy();
        return 0;
    }
    if (strcmp(name, "setters_return_ok") == 0) {
        test_setters_return_ok();
        return 0;
    }
    if (strcmp(name, "setters_reject_null_context") == 0) {
        test_setters_reject_null_context();
        return 0;
    }
    return 1;
}

int main(int argc, char *argv[]) {
    if (argc == 2)
        return run_named_test(argv[1]);

    test_version_string_defined();
    test_version_components_non_negative();
    test_version_components_in_string();
    test_context_create_destroy();
    test_setters_return_ok();
    test_setters_reject_null_context();
    return 0;
}
