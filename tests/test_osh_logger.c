#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "osh_logger.h"

void test_logger_basic() {
    const char *logfile = "test_logger.log";
    char line[512];
    FILE *f;
    int found_info = 0;
    int found_warn = 0;

    /* Setup logger */
    assert(osh_setup_logger(logfile, OSH_LOG_INFO) == 1);

    /* Test log level */
    assert(osh_get_loglevel() == OSH_LOG_INFO);
    assert(osh_set_loglevel(OSH_LOG_WARN) == 0);
    assert(osh_get_loglevel() == OSH_LOG_WARN);

    /* Write some logs */
    osh_info("Info message: %d", 42);
    osh_warn("Warning message: %s", "be careful");

    /* Close logger */
    assert(osh_close_logger() == 0);

    /* Check if file exists and contains something */
    f = fopen(logfile, "r");
    assert(f != NULL);


    while (fgets(line, sizeof(line), f)) {
        if (strstr(line, "Info message")) found_info = 1;
        if (strstr(line, "Warning message")) found_warn = 1;
    }
    fclose(f);

    /* Written logs alsway has info level irrespective of log level */
    assert(found_info == 1);
    assert(found_warn == 1);

    /* Clean up */
    remove(logfile);
}


int main() {
    test_logger_basic();
    printf("Logger test passed.\n");
    return 0;
}

