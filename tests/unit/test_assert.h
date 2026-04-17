#ifndef OSH_TEST_ASSERT_H
#define OSH_TEST_ASSERT_H

#include <stdio.h>
#include <stdlib.h>

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT_TRUE failed: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                            \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

#undef assert
#define assert(cond) ASSERT_TRUE(cond)

#endif /* OSH_TEST_ASSERT_H */
