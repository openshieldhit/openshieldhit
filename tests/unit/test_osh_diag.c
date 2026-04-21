#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "openshieldhit/diag.h"
#include "test_assert.h"

struct capture_diag {
    int calls;
    int last_level;
    char last_message[1024];
};

static void capture_emit(void *user, int level, char const *file, int line, char const *function, char const *message);
static void test_null_sink_is_noop(void);
static void test_null_emit_is_noop(void);
static void test_min_level_filters(void);
static void test_short_message_is_forwarded(void);
static void test_long_message_uses_heap_path(void);

int main(void) {
    test_null_sink_is_noop();
    test_null_emit_is_noop();
    test_min_level_filters();
    test_short_message_is_forwarded();
    test_long_message_uses_heap_path();
    return 0;
}

static void capture_emit(void *user, int level, char const *file, int line, char const *function, char const *message) {
    struct capture_diag *cap = (struct capture_diag *) user;
    (void) file;
    (void) line;
    (void) function;

    cap->calls += 1;
    cap->last_level = level;
    snprintf(cap->last_message, sizeof(cap->last_message), "%s", message ? message : "");
}

static void test_null_sink_is_noop(void) {
    osh_diag_emitf(NULL, OSH_DIAG_LEVEL_INFO, __FILE__, __LINE__, __func__, "hello %d", 1);
}

static void test_null_emit_is_noop(void) {
    struct osh_diag_sink diag = {0};

    diag.emit = NULL;
    diag.user = NULL;
    diag.min_level = OSH_DIAG_LEVEL_INFO;

    osh_diag_emitf(&diag, OSH_DIAG_LEVEL_INFO, __FILE__, __LINE__, __func__, "hello %d", 2);
}

static void test_min_level_filters(void) {
    struct capture_diag cap = {0};
    struct osh_diag_sink diag;

    diag.emit = capture_emit;
    diag.user = &cap;
    diag.min_level = OSH_DIAG_LEVEL_WARN;

    osh_diag_emitf(&diag, OSH_DIAG_LEVEL_INFO, __FILE__, __LINE__, __func__, "filtered");
    ASSERT_TRUE(cap.calls == 0);

    osh_diag_emitf(&diag, OSH_DIAG_LEVEL_ERROR, __FILE__, __LINE__, __func__, "visible");
    ASSERT_TRUE(cap.calls == 1);
    ASSERT_TRUE(cap.last_level == OSH_DIAG_LEVEL_ERROR);
    ASSERT_TRUE(strcmp(cap.last_message, "visible") == 0);
}

static void test_short_message_is_forwarded(void) {
    struct capture_diag cap = {0};
    struct osh_diag_sink diag;

    diag.emit = capture_emit;
    diag.user = &cap;
    diag.min_level = OSH_DIAG_LEVEL_INFO;

    osh_diag_emitf(&diag, OSH_DIAG_LEVEL_INFO, __FILE__, __LINE__, __func__, "beam %d", 42);

    ASSERT_TRUE(cap.calls == 1);
    ASSERT_TRUE(cap.last_level == OSH_DIAG_LEVEL_INFO);
    ASSERT_TRUE(strcmp(cap.last_message, "beam 42") == 0);
}

static void test_long_message_uses_heap_path(void) {
    struct capture_diag cap = {0};
    struct osh_diag_sink diag;
    char longbuf[900];
    size_t i;

    for (i = 0u; i < sizeof(longbuf) - 1u; ++i) {
        longbuf[i] = (char) ('a' + (i % 26u));
    }
    longbuf[sizeof(longbuf) - 1u] = '\0';

    diag.emit = capture_emit;
    diag.user = &cap;
    diag.min_level = OSH_DIAG_LEVEL_DEBUG;

    osh_diag_emitf(&diag, OSH_DIAG_LEVEL_DEBUG, __FILE__, __LINE__, __func__, "%s", longbuf);

    ASSERT_TRUE(cap.calls == 1);
    ASSERT_TRUE(cap.last_level == OSH_DIAG_LEVEL_DEBUG);
    ASSERT_TRUE(strcmp(cap.last_message, longbuf) == 0);
}
