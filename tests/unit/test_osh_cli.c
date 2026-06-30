#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cli/osh_cli.h"
#include "common/osh_duration.h"

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
        "openshieldhit", "-v", "-n", "1000", "--geo=/shared/geo.dat", "--beam=beam_override.dat", "run_17", NULL};

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

static void test_seedoffset_override(void) {
    char err[256];
    struct osh_cli_options opt;
    char *argv[] = {"openshieldhit", "-N", "123", "--seedoffset=456", NULL};

    int rc = osh_cli_parse(4, argv, &opt, err, sizeof(err));
    ASSERT_TRUE(rc == 0);
    ASSERT_TRUE(opt.has_seed_offset == 1);
    ASSERT_TRUE(opt.seed_offset == 456ULL);
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
    char *argv_minus[] = {"openshieldhit", "-n", "-1", NULL};

    int rc = osh_cli_parse(2, argv_plus, &opt, err, sizeof(err));
    ASSERT_TRUE(rc != 0);
    ASSERT_TRUE(strstr(err, "invalid integer value") != NULL);

    rc = osh_cli_parse(3, argv_minus, &opt, err, sizeof(err));
    ASSERT_TRUE(rc != 0);
    ASSERT_TRUE(strstr(err, "invalid integer value") != NULL);
}

static void test_seedoffset_rejects_sign_prefix(void) {
    char err[256];
    struct osh_cli_options opt;
    char *argv_plus[] = {"openshieldhit", "--seedoffset=+5", NULL};
    char *argv_minus[] = {"openshieldhit", "-N", "-1", NULL};

    int rc = osh_cli_parse(2, argv_plus, &opt, err, sizeof(err));
    ASSERT_TRUE(rc != 0);
    ASSERT_TRUE(strstr(err, "invalid integer value") != NULL);

    rc = osh_cli_parse(3, argv_minus, &opt, err, sizeof(err));
    ASSERT_TRUE(rc != 0);
    ASSERT_TRUE(strstr(err, "invalid integer value") != NULL);
}

static void test_seedoffset_rejects_above_9999(void) {
    char err[256];
    struct osh_cli_options opt;
    char *argv[] = {"openshieldhit", "--seedoffset=10000", NULL};

    int rc = osh_cli_parse(2, argv, &opt, err, sizeof(err));
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

static void test_profile_option(void) {
    char err[256];
    struct osh_cli_options opt;
    char *argv[] = {"openshieldhit", "--profile", "out/prof.json", "run_17", NULL};
    char *argv_eq[] = {"openshieldhit", "--profile=prof.json", NULL};
    char *argv_plain[] = {"openshieldhit", "run_17", NULL};
    char *argv_missing[] = {"openshieldhit", "--profile", NULL};

    int rc = osh_cli_parse(4, argv, &opt, err, sizeof(err));
    ASSERT_TRUE(rc == 0);
    ASSERT_TRUE(opt.profile_path != NULL && strcmp(opt.profile_path, "out/prof.json") == 0);
    ASSERT_TRUE(opt.workdir != NULL && strcmp(opt.workdir, "run_17") == 0);

    rc = osh_cli_parse(2, argv_eq, &opt, err, sizeof(err));
    ASSERT_TRUE(rc == 0);
    ASSERT_TRUE(opt.profile_path != NULL && strcmp(opt.profile_path, "prof.json") == 0);

    rc = osh_cli_parse(2, argv_plain, &opt, err, sizeof(err));
    ASSERT_TRUE(rc == 0);
    ASSERT_TRUE(opt.profile_path == NULL);

    rc = osh_cli_parse(2, argv_missing, &opt, err, sizeof(err));
    ASSERT_TRUE(rc != 0);
}

static void test_pool_capacity_option(void) {
    char err[256];
    struct osh_cli_options opt;
    char *argv[] = {"openshieldhit", "--pool-capacity", "4096", "run_17", NULL};
    char *argv_eq[] = {"openshieldhit", "--pool-capacity=256", NULL};
    char *argv_zero[] = {"openshieldhit", "--pool-capacity=0", NULL};
    char *argv_plain[] = {"openshieldhit", "run_17", NULL};
    char *argv_bad[] = {"openshieldhit", "--pool-capacity=abc", NULL};
    char *argv_missing[] = {"openshieldhit", "--pool-capacity", NULL};

    int rc = osh_cli_parse(4, argv, &opt, err, sizeof(err));
    ASSERT_TRUE(rc == 0);
    ASSERT_TRUE(opt.has_pool_capacity == 1);
    ASSERT_TRUE(opt.pool_capacity == 4096ULL);
    ASSERT_TRUE(opt.workdir != NULL && strcmp(opt.workdir, "run_17") == 0);

    rc = osh_cli_parse(2, argv_eq, &opt, err, sizeof(err));
    ASSERT_TRUE(rc == 0);
    ASSERT_TRUE(opt.has_pool_capacity == 1);
    ASSERT_TRUE(opt.pool_capacity == 256ULL);

    /* 0 is valid and means "use the compiled default" — must not be rejected. */
    rc = osh_cli_parse(2, argv_zero, &opt, err, sizeof(err));
    ASSERT_TRUE(rc == 0);
    ASSERT_TRUE(opt.has_pool_capacity == 1);
    ASSERT_TRUE(opt.pool_capacity == 0ULL);

    /* Absent flag leaves the override off. */
    rc = osh_cli_parse(2, argv_plain, &opt, err, sizeof(err));
    ASSERT_TRUE(rc == 0);
    ASSERT_TRUE(opt.has_pool_capacity == 0);

    rc = osh_cli_parse(2, argv_bad, &opt, err, sizeof(err));
    ASSERT_TRUE(rc != 0);
    ASSERT_TRUE(strstr(err, "invalid integer value") != NULL);

    rc = osh_cli_parse(2, argv_missing, &opt, err, sizeof(err));
    ASSERT_TRUE(rc != 0);
}

static void test_parse_duration_good(void) {
    double s = -1.0;

    /* Bare number = seconds. */
    ASSERT_TRUE(osh_parse_duration("0", &s) == 1 && s == 0.0);
    ASSERT_TRUE(osh_parse_duration("500", &s) == 1 && s == 500.0);

    /* Unit suffixes (case-insensitive). */
    ASSERT_TRUE(osh_parse_duration("30s", &s) == 1 && s == 30.0);
    ASSERT_TRUE(osh_parse_duration("30S", &s) == 1 && s == 30.0);
    ASSERT_TRUE(osh_parse_duration("30m", &s) == 1 && s == 1800.0);
    ASSERT_TRUE(osh_parse_duration("1h", &s) == 1 && s == 3600.0);
    ASSERT_TRUE(osh_parse_duration("1H", &s) == 1 && s == 3600.0);

    /* Fractional values are allowed. */
    ASSERT_TRUE(osh_parse_duration("1.5h", &s) == 1 && s == 5400.0);
    ASSERT_TRUE(osh_parse_duration("0.5", &s) == 1 && s == 0.5);
}

static void test_parse_duration_bad(void) {
    double s = -1.0;

    ASSERT_TRUE(osh_parse_duration(NULL, &s) == 0);
    ASSERT_TRUE(osh_parse_duration("", &s) == 0);
    ASSERT_TRUE(osh_parse_duration("abc", &s) == 0);
    ASSERT_TRUE(osh_parse_duration("30x", &s) == 0);  /* unknown unit */
    ASSERT_TRUE(osh_parse_duration("30ss", &s) == 0); /* trailing junk after unit */
    ASSERT_TRUE(osh_parse_duration(" 30", &s) == 0);  /* leading whitespace rejected */
    ASSERT_TRUE(osh_parse_duration("+30", &s) == 0);  /* sign prefix rejected */
    ASSERT_TRUE(osh_parse_duration("-30", &s) == 0);  /* negative rejected */
    ASSERT_TRUE(osh_parse_duration("30 ", &s) == 0);  /* trailing whitespace rejected */
}

static void test_max_time_option(void) {
    char err[256];
    struct osh_cli_options opt;
    char *argv[] = {"openshieldhit", "--max-time", "30m", "run_17", NULL};
    char *argv_eq[] = {"openshieldhit", "--max-time=1h", NULL};
    char *argv_plain[] = {"openshieldhit", "run_17", NULL};
    char *argv_bad[] = {"openshieldhit", "--max-time=nope", NULL};
    char *argv_missing[] = {"openshieldhit", "--max-time", NULL};

    int rc = osh_cli_parse(4, argv, &opt, err, sizeof(err));
    ASSERT_TRUE(rc == 0);
    ASSERT_TRUE(opt.has_max_time == 1);
    ASSERT_TRUE(opt.max_time_s == 1800.0);
    ASSERT_TRUE(opt.workdir != NULL && strcmp(opt.workdir, "run_17") == 0);

    rc = osh_cli_parse(2, argv_eq, &opt, err, sizeof(err));
    ASSERT_TRUE(rc == 0);
    ASSERT_TRUE(opt.has_max_time == 1);
    ASSERT_TRUE(opt.max_time_s == 3600.0);

    /* Absent flag leaves the override off. */
    rc = osh_cli_parse(2, argv_plain, &opt, err, sizeof(err));
    ASSERT_TRUE(rc == 0);
    ASSERT_TRUE(opt.has_max_time == 0);

    rc = osh_cli_parse(2, argv_bad, &opt, err, sizeof(err));
    ASSERT_TRUE(rc != 0);
    ASSERT_TRUE(strstr(err, "invalid duration value") != NULL);

    rc = osh_cli_parse(2, argv_missing, &opt, err, sizeof(err));
    ASSERT_TRUE(rc != 0);
}

static int run_named_test(char const *name) {
    if (strcmp(name, "parse_duration_good") == 0) {
        test_parse_duration_good();
        return 0;
    }
    if (strcmp(name, "parse_duration_bad") == 0) {
        test_parse_duration_bad();
        return 0;
    }
    if (strcmp(name, "max_time_option") == 0) {
        test_max_time_option();
        return 0;
    }
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
    if (strcmp(name, "seedoffset_override") == 0) {
        test_seedoffset_override();
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
    if (strcmp(name, "seedoffset_rejects_sign_prefix") == 0) {
        test_seedoffset_rejects_sign_prefix();
        return 0;
    }
    if (strcmp(name, "seedoffset_rejects_above_9999") == 0) {
        test_seedoffset_rejects_above_9999();
        return 0;
    }
    if (strcmp(name, "rejects_extra_positional_argument") == 0) {
        test_rejects_extra_positional_argument();
        return 0;
    }
    if (strcmp(name, "profile_option") == 0) {
        test_profile_option();
        return 0;
    }
    if (strcmp(name, "pool_capacity_option") == 0) {
        test_pool_capacity_option();
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
    test_seedoffset_override();
    test_nstat_rejects_leading_whitespace();
    test_nstat_rejects_sign_prefix();
    test_seedoffset_rejects_sign_prefix();
    test_seedoffset_rejects_above_9999();
    test_rejects_extra_positional_argument();
    test_profile_option();
    test_pool_capacity_option();
    test_parse_duration_good();
    test_parse_duration_bad();
    test_max_time_option();
    return 0;
}
