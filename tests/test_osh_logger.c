#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_logger.h"

static int file_contains(const char *path, const char *needle) {
    char line[1024];
    FILE *f = fopen(path, "r");
    if (!f)
        return 0;

    while (fgets(line, sizeof(line), f)) {
        if (strstr(line, needle)) {
            fclose(f);
            return 1;
        }
    }
    fclose(f);
    return 0;
}

static int file_not_contains(const char *path, const char *needle) {
    return !file_contains(path, needle);
}

void test_logger_basic(void) {
    const char *logfile = "test_logger.log";

    /* Start clean */
    remove(logfile);

    /* Init default logger */
    assert(osh_log_init(OSH_LOG_INFO, OSH_LOG_F_TIMESTAMP | OSH_LOG_F_FILELINE) == 0);

    /* Add file sink */
    assert(osh_log_add_file(logfile, /*append=*/0) == 0);

    /* Level getters/setters */
    assert(osh_log_get_level() == OSH_LOG_INFO);
    assert(osh_log_set_level(OSH_LOG_WARN) == 0);
    assert(osh_log_get_level() == OSH_LOG_WARN);

    /* With level WARN: INFO should be suppressed, WARN should appear */
    osh_info("Info message (should NOT appear): %d", 42);
    osh_warn("Warning message A: %s", "be careful");

    /* Lower threshold to INFO so INFO appears */
    assert(osh_log_set_level(OSH_LOG_INFO) == 0);
    osh_info("Info message (should appear): %d", 43);
    osh_warn("Warning message B: %s", "still careful");

    /* Ensure everything is written */
    osh_log_flush();
    osh_log_close();

    /* Verify file exists */
    {
        FILE *f = fopen(logfile, "r");
        assert(f != NULL);
        fclose(f);
    }

    /* Verify contents */
    assert(file_contains(logfile, "Warning message A"));
    assert(file_contains(logfile, "Warning message B"));
    assert(file_contains(logfile, "Info message (should appear)"));
    assert(file_not_contains(logfile, "Info message (should NOT appear)"));

    /* Clean up */
    remove(logfile);
}

int main(void) {
    test_logger_basic();
    printf("Logger test passed.\n");
    return 0;
}
