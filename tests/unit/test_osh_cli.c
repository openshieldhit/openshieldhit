#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cli/osh_cli.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

static void test_version_short_flag(void) {
    char err[256];
    struct osh_cli_options opt;
    char *argv[] = {"openshieldhit", "-V", NULL};

    int rc = osh_cli_parse(2, argv, &opt, err, sizeof(err));
    ASSERT_TRUE(rc == 0);
    ASSERT_TRUE(opt.action == OSH_CLI_ACTION_VERSION);
    ASSERT_TRUE(opt.verbose == 0);
}

static void test_verbose_and_dir_options(void) {
    char err[256];
    struct osh_cli_options opt;
    char *argv[] = {"openshieldhit", "-v", "--dry-run", "--workdir", "runs/case01", "--outdir", "/tmp/out", NULL};

    int rc = osh_cli_parse(7, argv, &opt, err, sizeof(err));
    ASSERT_TRUE(rc == 0);
    ASSERT_TRUE(opt.action == OSH_CLI_ACTION_RUN);
    ASSERT_TRUE(opt.dry_run == 1);
    ASSERT_TRUE(opt.verbose == 1);
    ASSERT_TRUE(opt.workdir != NULL && strcmp(opt.workdir, "runs/case01") == 0);
    ASSERT_TRUE(opt.out_dir != NULL && strcmp(opt.out_dir, "/tmp/out") == 0);
}

static void test_nstat_and_overrides(void) {
    char err[256];
    struct osh_cli_options opt;
    char *argv[] = {
        "openshieldhit", "-v", "-N", "1000", "--geo=/shared/geo.dat", "--beam=beam_override.dat", "run_17", NULL};

    int rc = osh_cli_parse(7, argv, &opt, err, sizeof(err));
    ASSERT_TRUE(rc == 0);
    ASSERT_TRUE(opt.action == OSH_CLI_ACTION_RUN);
    ASSERT_TRUE(opt.verbose == 1);
    ASSERT_TRUE(opt.has_nstat == 1);
    ASSERT_TRUE(opt.nstat == 1000ULL);
    ASSERT_TRUE(opt.workdir != NULL && strcmp(opt.workdir, "run_17") == 0);
    ASSERT_TRUE(opt.geo_path != NULL && strcmp(opt.geo_path, "/shared/geo.dat") == 0);
    ASSERT_TRUE(opt.beam_path != NULL && strcmp(opt.beam_path, "beam_override.dat") == 0);
}

static void test_nstat_rejects_leading_whitespace(void) {
    char err[256];
    struct osh_cli_options opt;
    char *argv[] = {"openshieldhit", "--nstat", " 123", NULL};

    int rc = osh_cli_parse(3, argv, &opt, err, sizeof(err));
    ASSERT_TRUE(rc != 0);
    ASSERT_TRUE(strstr(err, "invalid integer value") != NULL);
}

static void test_nstat_rejects_sign_prefix(void) {
    char err[256];
    struct osh_cli_options opt;
    char *argv_plus[] = {"openshieldhit", "--nstat=+5", NULL};
    char *argv_minus[] = {"openshieldhit", "-N", "-1", NULL};

    int rc = osh_cli_parse(2, argv_plus, &opt, err, sizeof(err));
    ASSERT_TRUE(rc != 0);
    ASSERT_TRUE(strstr(err, "invalid integer value") != NULL);

    rc = osh_cli_parse(3, argv_minus, &opt, err, sizeof(err));
    ASSERT_TRUE(rc != 0);
    ASSERT_TRUE(strstr(err, "invalid integer value") != NULL);
}

static void test_rejects_extra_positional_argument(void) {
    char err[256];
    struct osh_cli_options opt;
    char *argv[] = {"openshieldhit", "run_a", "run_b", NULL};

    int rc = osh_cli_parse(3, argv, &opt, err, sizeof(err));
    ASSERT_TRUE(rc != 0);
    ASSERT_TRUE(strstr(err, "unexpected positional argument") != NULL);
}

static int run_named_test(char const *name) {
    if (strcmp(name, "version_short_flag") == 0) {
        test_version_short_flag();
        return 0;
    }
    if (strcmp(name, "verbose_and_dir_options") == 0) {
        test_verbose_and_dir_options();
        return 0;
    }
    if (strcmp(name, "nstat_and_overrides") == 0) {
        test_nstat_and_overrides();
        return 0;
    }
    if (strcmp(name, "nstat_rejects_leading_whitespace") == 0) {
        test_nstat_rejects_leading_whitespace();
        return 0;
    }
    if (strcmp(name, "nstat_rejects_sign_prefix") == 0) {
        test_nstat_rejects_sign_prefix();
        return 0;
    }
    if (strcmp(name, "rejects_extra_positional_argument") == 0) {
        test_rejects_extra_positional_argument();
        return 0;
    }
    return 1;
}

int main(int argc, char *argv[]) {
    if (argc == 2) {
        return run_named_test(argv[1]);
    }

    test_version_short_flag();
    test_verbose_and_dir_options();
    test_nstat_and_overrides();
    test_nstat_rejects_leading_whitespace();
    test_nstat_rejects_sign_prefix();
    test_rejects_extra_positional_argument();
    return 0;
}
