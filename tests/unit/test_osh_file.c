#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_file.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

static void test_relative_path_join(void);
static void test_dirname_accepts_backslashes(void);
static void test_windows_absolute_path_detection(void);
static void test_remove_dir_removes_files_and_directory(void);
static void test_remove_dir_rejects_invalid_input(void);

static int _count_entry(char const *path, void *user);

int main(void) {
    test_relative_path_join();
    test_dirname_accepts_backslashes();
    test_windows_absolute_path_detection();
    test_remove_dir_removes_files_and_directory();
    test_remove_dir_rejects_invalid_input();
    return 0;
}

static void test_relative_path_join(void) {
    char *resolved = NULL;

    ASSERT_TRUE(osh_relative_path_to_file(&resolved, "/tmp/work", "out.dat") == 0);
    ASSERT_TRUE(strcmp(resolved, "/tmp/work/out.dat") == 0);
    free(resolved);
}

static void test_dirname_accepts_backslashes(void) {
    char *dir;

    dir = osh_path_dirname("C:\\work\\detect.dat");
    ASSERT_TRUE(dir != NULL);
    ASSERT_TRUE(strcmp(dir, "C:\\work") == 0);
    free(dir);
}

static void test_windows_absolute_path_detection(void) {
    char *resolved = NULL;

#if defined(_WIN32)
    ASSERT_TRUE(osh_relative_path_to_file(&resolved, "C:/base", "C:/tmp/out.dat") == 0);
    ASSERT_TRUE(strcmp(resolved, "C:/tmp/out.dat") == 0);
    free(resolved);
    resolved = NULL;

    ASSERT_TRUE(osh_relative_path_to_file(&resolved, "C:/base", "C:\\tmp\\out.dat") == 0);
    ASSERT_TRUE(strcmp(resolved, "C:\\tmp\\out.dat") == 0);
    free(resolved);
    resolved = NULL;

    ASSERT_TRUE(osh_relative_path_to_file(&resolved, "C:/base", "C:tmp\\out.dat") == 0);
    ASSERT_TRUE(strcmp(resolved, "C:/base/C:tmp\\out.dat") == 0);
    free(resolved);
#else
    ASSERT_TRUE(osh_relative_path_to_file(&resolved, "/tmp/base", "C:/tmp/out.dat") == 0);
    ASSERT_TRUE(strcmp(resolved, "/tmp/base/C:/tmp/out.dat") == 0);
    free(resolved);
#endif
}

static int _count_entry(char const *path, void *user) {
    int *count = (int *) user;

    (void) path;
    (*count)++;
    return 1;
}

/* osh_path_remove_dir() is the teardown counterpart to osh_path_ensure_dir(),
 * used by tests (e.g. test_osh_run_dump.c) to avoid leaving scratch output
 * directories behind. Cover both the create-populate-remove path and the
 * missing-path no-op. */
static void test_remove_dir_removes_files_and_directory(void) {
    char const *dir = "test_osh_file_remove_dir_scratch";
    char file_path[256];
    FILE *fp;
    int count;

    ASSERT_TRUE(osh_path_ensure_dir(dir) == OSH_OK);
    snprintf(file_path, sizeof(file_path), "%s/scratch.tmp", dir);
    fp = fopen(file_path, "w");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fputs("x", fp) >= 0);
    ASSERT_TRUE(fclose(fp) == 0);

    ASSERT_TRUE(osh_path_remove_dir(dir) == OSH_OK);

    count = 0;
    ASSERT_TRUE(osh_dir_foreach_file(dir, _count_entry, &count) != OSH_OK);
    ASSERT_TRUE(count == 0);

    /* Already gone: a repeat call is a no-op, not an error. */
    ASSERT_TRUE(osh_path_remove_dir(dir) == OSH_OK);
}

/* Invalid arguments and a non-directory path are rejected rather than
 * silently treated as "nothing to do". */
static void test_remove_dir_rejects_invalid_input(void) {
    char const *file_path = "test_osh_file_remove_dir_not_a_dir.tmp";
    FILE *fp;

    ASSERT_TRUE(osh_path_remove_dir(NULL) == OSH_EINVAL);
    ASSERT_TRUE(osh_path_remove_dir("") == OSH_EINVAL);

    fp = fopen(file_path, "w");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fputs("x", fp) >= 0);
    ASSERT_TRUE(fclose(fp) == 0);

    ASSERT_TRUE(osh_path_remove_dir(file_path) == OSH_EIO);

    remove(file_path);
}
