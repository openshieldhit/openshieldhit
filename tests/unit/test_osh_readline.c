#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_file.h"
#include "common/osh_readline.h"

#ifndef OSH_TEST_FIXTURES_DIR
#error "OSH_TEST_FIXTURES_DIR must be defined by the build system"
#endif
#define TEST_PATH OSH_TEST_FIXTURES_DIR "/test01/beam.dat"

void test_readline_key(char const *filename) {
    struct oshfile *oshf;

    char *lline = NULL;
    char *key = NULL;
    char *args = NULL;
    int lineno;
    int len;

    oshf = osh_fopen(filename);
    assert(oshf != NULL);

    len = osh_readline_key(oshf, &lline, &key, &args, &lineno);
    assert(len > 0);
    assert(key != NULL);
    assert(args != NULL);

    assert(strcmp(key, "RNDSEED") == 0);
    assert(strcmp(args, "89736501") == 0);
    assert(lineno == 1);

    free(lline);
    osh_fclose(oshf);

    printf("test_readline_key passed.\n");
}

void test_count_keys(char const *filename) {

    struct oshfile *oshf;
    int count = 0;
    char *line = NULL;
    char *key = NULL;
    char *args = NULL;
    int lineno;

    oshf = osh_fopen(filename);
    assert(oshf != NULL);

    while (osh_readline_key(oshf, &line, &key, &args, &lineno) != -1) {
        assert(key != NULL); /* Every valid line must have a key   */
        count++;
        assert(args != NULL); /* Every valid line must have args    */
        assert(lineno > 0);   /* Every valid line must have a lineno */
        assert(key != NULL);  /* Every valid key must be non-NULL */
/* print the key and the argument */
#ifdef DEBUG
        printf("Found key in line %d: %s with args: %s\n", lineno, key, args);
#endif
        free(line);
    }

    osh_fclose(oshf);
    assert(count == 13);
    printf("test_count_keys() passed. Found %d keys.\n", count);
}

int main(void) {
    test_readline_key(TEST_PATH);
    test_count_keys(TEST_PATH);
    return 0;
}
