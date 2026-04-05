#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "openshieldhit/openshieldhit.h"

/* Set by main() when a verbosity level is passed as the first argument.
 * Tests that call openshieldhit_run() propagate this to cfg.log_level so the
 * library initialises the logger at the requested level. */
static int g_verbosity = 0;

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

static FILE *osh_test_tmpfile(void) {
#if defined(_MSC_VER)
    FILE *fp = NULL;
    if (tmpfile_s(&fp) != 0) {
        return NULL;
    }
    return fp;
#else
    return tmpfile();
#endif
}

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

static void test_config_init_sets_size(void) {
    openshieldhit_config_t cfg = OPENSHIELDHIT_CONFIG_INIT;

    ASSERT_TRUE(cfg.size == sizeof(openshieldhit_config_t));
    ASSERT_TRUE(cfg.workdir == NULL);
    ASSERT_TRUE(cfg.run_mode == OPENSHIELDHIT_RUN_NORMAL);
    ASSERT_TRUE(cfg.has_nstat == 0);
}

static void test_configure_returns_ok(void) {
    openshieldhit_context_t *ctx = openshieldhit_context_create();
    openshieldhit_config_t cfg = OPENSHIELDHIT_CONFIG_INIT;

    ASSERT_TRUE(ctx != NULL);

    cfg.workdir = "/tmp";
    cfg.out_dir = "/tmp/out";
    cfg.geo_path = "geo.dat";
    cfg.beam_path = "beam.dat";
    cfg.mat_path = "mat.dat";
    cfg.detect_path = "detect.dat";
    cfg.nstat = 1000ULL;
    cfg.has_nstat = 1;
    cfg.log_level = g_verbosity;
    cfg.run_mode = OPENSHIELDHIT_RUN_VALIDATE;

    ASSERT_TRUE(openshieldhit_context_configure(ctx, &cfg) == OPENSHIELDHIT_STATUS_OK);

    openshieldhit_context_destroy(ctx);
}

static void test_configure_rejects_invalid_arguments(void) {
    openshieldhit_context_t *ctx;
    openshieldhit_config_t cfg = OPENSHIELDHIT_CONFIG_INIT;

    ASSERT_TRUE(openshieldhit_context_configure(NULL, &cfg) == OPENSHIELDHIT_STATUS_INVALID_ARGUMENT);

    ctx = openshieldhit_context_create();
    ASSERT_TRUE(ctx != NULL);
    cfg.size = 0;
    ASSERT_TRUE(openshieldhit_context_configure(ctx, &cfg) == OPENSHIELDHIT_STATUS_INVALID_ARGUMENT);
    openshieldhit_context_destroy(ctx);
}

static void test_last_error_empty_on_fresh_context(void) {
    openshieldhit_context_t *ctx = openshieldhit_context_create();
    ASSERT_TRUE(ctx != NULL);
    ASSERT_TRUE(openshieldhit_last_error(ctx) != NULL);
    ASSERT_TRUE(openshieldhit_last_error(ctx)[0] == '\0');
    openshieldhit_context_destroy(ctx);
}

static void test_last_error_empty_on_null_context(void) {
    char const *err = openshieldhit_last_error(NULL);
    ASSERT_TRUE(err != NULL);
    ASSERT_TRUE(err[0] == '\0');
}

static void test_last_error_set_on_unsupported_run_mode(void) {
    openshieldhit_context_t *ctx = openshieldhit_context_create();
    openshieldhit_config_t cfg = OPENSHIELDHIT_CONFIG_INIT;
    enum openshieldhit_status rc;

    ASSERT_TRUE(ctx != NULL);
    /* NORMAL mode is not yet implemented — run() should fail and set last_error. */
    cfg.run_mode = OPENSHIELDHIT_RUN_NORMAL;
    ASSERT_TRUE(openshieldhit_context_configure(ctx, &cfg) == OPENSHIELDHIT_STATUS_OK);

    rc = openshieldhit_run(ctx, NULL, NULL);
    ASSERT_TRUE(rc == OPENSHIELDHIT_STATUS_NOT_SUPPORTED);
    ASSERT_TRUE(openshieldhit_last_error(ctx)[0] != '\0');

    openshieldhit_context_destroy(ctx);
}

static void test_validate_applies_nstat_override(void) {
    openshieldhit_context_t *ctx = openshieldhit_context_create();
    openshieldhit_config_t cfg = OPENSHIELDHIT_CONFIG_INIT;
    enum openshieldhit_status rc;
    FILE *out;
    char buf[4096];
    char workdir[512];
    size_t nread;

    ASSERT_TRUE(ctx != NULL);

    snprintf(workdir, sizeof(workdir), "%s/tests/cases/02_sobp", OSH_PROJECT_SOURCE_DIR);
    cfg.workdir = workdir;
    cfg.run_mode = OPENSHIELDHIT_RUN_VALIDATE;
    cfg.nstat = 1234ULL;
    cfg.has_nstat = 1;
    cfg.log_level = g_verbosity;

    ASSERT_TRUE(openshieldhit_context_configure(ctx, &cfg) == OPENSHIELDHIT_STATUS_OK);

    out = osh_test_tmpfile();
    ASSERT_TRUE(out != NULL);

    rc = openshieldhit_run(ctx, out, NULL);
    ASSERT_TRUE(rc == OPENSHIELDHIT_STATUS_OK);

    ASSERT_TRUE(fseek(out, 0, SEEK_SET) == 0);
    nread = fread(buf, 1, sizeof(buf) - 1u, out);
    buf[nread] = '\0';

    ASSERT_TRUE(strstr(buf, "Requested nstat  : 1234") != NULL);
    ASSERT_TRUE(strstr(buf, "Applied nstat override: 1234") != NULL);

    ASSERT_TRUE(fclose(out) == 0);
    openshieldhit_context_destroy(ctx);
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
    if (strcmp(name, "config_init_sets_size") == 0) {
        test_config_init_sets_size();
        return 0;
    }
    if (strcmp(name, "configure_returns_ok") == 0) {
        test_configure_returns_ok();
        return 0;
    }
    if (strcmp(name, "configure_rejects_invalid_arguments") == 0) {
        test_configure_rejects_invalid_arguments();
        return 0;
    }
    if (strcmp(name, "last_error_empty_on_fresh_context") == 0) {
        test_last_error_empty_on_fresh_context();
        return 0;
    }
    if (strcmp(name, "last_error_empty_on_null_context") == 0) {
        test_last_error_empty_on_null_context();
        return 0;
    }
    if (strcmp(name, "last_error_set_on_unsupported_run_mode") == 0) {
        test_last_error_set_on_unsupported_run_mode();
        return 0;
    }
    if (strcmp(name, "validate_applies_nstat_override") == 0) {
        test_validate_applies_nstat_override();
        return 0;
    }
    return 1;
}

int main(int argc, char *argv[]) {
    if (argc == 3) {
        g_verbosity = atoi(argv[1]);
        return run_named_test(argv[2]);
    }
    if (argc == 2) {
        return run_named_test(argv[1]);
    }

    test_version_string_defined();
    test_version_components_non_negative();
    test_version_components_in_string();
    test_context_create_destroy();
    test_config_init_sets_size();
    test_configure_returns_ok();
    test_configure_rejects_invalid_arguments();
    test_last_error_empty_on_fresh_context();
    test_last_error_empty_on_null_context();
    test_last_error_set_on_unsupported_run_mode();
    test_validate_applies_nstat_override();
    return 0;
}
