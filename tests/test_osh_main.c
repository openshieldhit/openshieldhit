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
    snprintf(
        expected,
        sizeof(expected),
        "%d.%d.%d",
        openshieldhit_version_major(),
        openshieldhit_version_minor(),
        openshieldhit_version_patch());
    ASSERT_TRUE(strstr(openshieldhit_version_string(), expected) != NULL);
}

static void test_cli_version_short_flag(void) {
    int rc;
    char err[256];
    struct openshieldhit_cli_options opt;
    char *argv[] = {"openshieldhit", "-V", NULL};

    rc = openshieldhit_cli_parse(2, argv, &opt, err, sizeof(err));
    ASSERT_TRUE(rc == 0);
    ASSERT_TRUE(opt.action == OPENSHIELDHIT_CLI_ACTION_VERSION);
    ASSERT_TRUE(opt.verbose == 0);
}

static void test_cli_verbose_and_dir_options(void) {
    int rc;
    char err[256];
    struct openshieldhit_cli_options opt;
    char *argv[] = {
        "openshieldhit", "-v", "--dry-run", "--workdir", "runs/case01", "--outdir", "/tmp/out", NULL};

    rc = openshieldhit_cli_parse(7, argv, &opt, err, sizeof(err));
    ASSERT_TRUE(rc == 0);
    ASSERT_TRUE(opt.action == OPENSHIELDHIT_CLI_ACTION_RUN);
    ASSERT_TRUE(opt.dry_run == 1);
    ASSERT_TRUE(opt.verbose == 1);
    ASSERT_TRUE(opt.workdir != NULL);
    ASSERT_TRUE(strcmp(opt.workdir, "runs/case01") == 0);
    ASSERT_TRUE(opt.out_dir != NULL);
    ASSERT_TRUE(strcmp(opt.out_dir, "/tmp/out") == 0);
}

static void test_cli_nstat_and_overrides(void) {
    int rc;
    char err[256];
    struct openshieldhit_cli_options opt;
    char *argv[] = {
        "openshieldhit",
        "-v",
        "-N",
        "1000",
        "--geo=/shared/geo.dat",
        "--beam=beam_override.dat",
        "run_17",
        NULL};

    rc = openshieldhit_cli_parse(7, argv, &opt, err, sizeof(err));
    ASSERT_TRUE(rc == 0);
    ASSERT_TRUE(opt.action == OPENSHIELDHIT_CLI_ACTION_RUN);
    ASSERT_TRUE(opt.verbose == 1);
    ASSERT_TRUE(opt.has_nstat == 1);
    ASSERT_TRUE(opt.nstat == 1000ULL);
    ASSERT_TRUE(opt.workdir != NULL);
    ASSERT_TRUE(strcmp(opt.workdir, "run_17") == 0);
    ASSERT_TRUE(opt.geo_path != NULL);
    ASSERT_TRUE(strcmp(opt.geo_path, "/shared/geo.dat") == 0);
    ASSERT_TRUE(opt.beam_path != NULL);
    ASSERT_TRUE(strcmp(opt.beam_path, "beam_override.dat") == 0);
}

static void test_cli_rejects_extra_positional_argument(void) {
    int rc;
    char err[256];
    struct openshieldhit_cli_options opt;
    char *argv[] = {"openshieldhit", "run_a", "run_b", NULL};

    rc = openshieldhit_cli_parse(3, argv, &opt, err, sizeof(err));
    ASSERT_TRUE(rc != 0);
    ASSERT_TRUE(strstr(err, "unexpected positional argument") != NULL);
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
    if (strcmp(name, "cli_version_short_flag") == 0) {
        test_cli_version_short_flag();
        return 0;
    }
    if (strcmp(name, "cli_verbose_and_dir_options") == 0) {
        test_cli_verbose_and_dir_options();
        return 0;
    }
    if (strcmp(name, "cli_nstat_and_overrides") == 0) {
        test_cli_nstat_and_overrides();
        return 0;
    }
    if (strcmp(name, "cli_rejects_extra_positional_argument") == 0) {
        test_cli_rejects_extra_positional_argument();
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
    test_cli_version_short_flag();
    test_cli_verbose_and_dir_options();
    test_cli_nstat_and_overrides();
    test_cli_rejects_extra_positional_argument();
    return 0;
}
