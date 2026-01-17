#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "osh_logger.h"

#if defined(_WIN32)
#include <windows.h>
#else
#include <pthread.h>
#include <unistd.h>
#endif

#ifndef THREADS
#define THREADS 4
#endif

#ifndef ITERS
#define ITERS 2000
#endif

struct thread_arg {
    int tid;
    int iters;
};

/* --- line counting helper --- */
static long count_newlines_in_file(const char *path) {
    FILE *f = fopen(path, "rb");
    long count = 0;
    int c;

    if (!f)
        return -1;

    while ((c = fgetc(f)) != EOF) {
        if (c == '\n')
            count++;
    }
    fclose(f);
    return count;
}

/* --- thread entry --- */
#if defined(_WIN32)

static DWORD WINAPI thread_main(LPVOID param) {
    struct thread_arg *a = (struct thread_arg *)param;
    int tid = a->tid;
    int iters = a->iters;

    for (int i = 0; i < iters; i++) {
        /* Include a predictable token to help manual inspection if needed */
        osh_logger_log_ex(osh_log_default(), OSH_LOG_INFO, 0u, __FILE__, __LINE__, __func__,
                          "T=%d i=%d token=THREADTEST", tid, i);
    }
    return 0;
}

#else

static void *thread_main(void *param) {
    struct thread_arg *a = (struct thread_arg *)param;
    int tid = a->tid;
    int iters = a->iters;

    for (int i = 0; i < iters; i++) {
        osh_logger_log_ex(osh_log_default(), OSH_LOG_INFO, 0u, __FILE__, __LINE__, __func__,
                          "T=%d i=%d token=THREADTEST", tid, i);
    }
    return NULL;
}

#endif

static void test_logger_threads(void) {
    const char *logfile = "test_logger_threads.log";
    const int nthreads = THREADS;
    const int iters = ITERS;

    remove(logfile);

    /* Init + file sink; keep console quiet */
    assert(osh_log_init(OSH_LOG_INFO, OSH_LOG_F_NONE) == 0);
    assert(osh_log_add_file(logfile, /*append=*/0) == 0);
    (void)osh_log_enable_stdout(0);

#if defined(_WIN32)
    HANDLE th[THREADS];
    DWORD tid[THREADS];
    struct thread_arg args[THREADS];

    for (int i = 0; i < nthreads; i++) {
        args[i].tid = i;
        args[i].iters = iters;
        th[i] = CreateThread(NULL, 0, thread_main, &args[i], 0, &tid[i]);
        assert(th[i] != NULL);
    }

    /* Wait for all */
    {
        DWORD rc = WaitForMultipleObjects((DWORD)nthreads, th, TRUE, INFINITE);
        assert(rc >= WAIT_OBJECT_0 && rc < WAIT_OBJECT_0 + (DWORD)nthreads);
    }

    for (int i = 0; i < nthreads; i++) {
        CloseHandle(th[i]);
    }

#else
    pthread_t th[THREADS];
    struct thread_arg args[THREADS];

    for (int i = 0; i < nthreads; i++) {
        args[i].tid = i;
        args[i].iters = iters;
        int rc = pthread_create(&th[i], NULL, thread_main, &args[i]);
        assert(rc == 0);
    }

    for (int i = 0; i < nthreads; i++) {
        int rc = pthread_join(th[i], NULL);
        assert(rc == 0);
    }
#endif

    osh_log_flush();
    osh_log_close();

    /* Verify file exists */
    {
        FILE *f = fopen(logfile, "rb");
        assert(f != NULL);
        fclose(f);
    }

    /* Verify exactly N*ITERS lines */
    {
        long nl = count_newlines_in_file(logfile);
        assert(nl == (long)nthreads * (long)iters);
    }

    //
}

int main(void) {
    test_logger_threads();
    printf("Logger threaded test passed. (%d threads x %d lines)\n", THREADS, ITERS);
    return 0;
}