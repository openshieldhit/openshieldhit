#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "openshieldhit/file.h"

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

int main(void) {
    test_relative_path_join();
    test_dirname_accepts_backslashes();
    test_windows_absolute_path_detection();
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
