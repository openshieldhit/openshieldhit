#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>

#include "osh_exit.h"
#include "osh_logger.h"

#define _MAX_LINE_LENGTH 4096   /* Consider making this global somehow */

/* help GCC to tell that this function takes printf-style arguments */
#if defined(__GNUC__) || defined(__clang__)
__attribute__((format(printf, 2, 3)))
__attribute__((noreturn))
#endif

static int _log_level = 4;          /* this way, there will be printouts, even if the logger was not setup first */
static char *_log_filename;
static int _logger_save_active = 0; /* init FALSE. If logger was setup with an output file, then enable this. */
static FILE *_fplog = NULL;

/* warning message database only visbile inside this .c file */
static struct message_db {
    char **msg; /* list of unique messages */
    int *count; /* number of times a message has been occuring */
    int *level; /* message level */
    int len;    /* number of unique messages */
} message_db;

/* internal functions of the warning message database */
void _message_db_init(void);
void _message_db_summary(void);
void _message_db_free(void);
int _message_db_check(char const *s, int level);


int osh_setup_logger(char const *path, int log_level) {

    int len;

    if (_logger_save_active) {
        osh_warn("osh_setup_logger() called while logger already active. Ignoring.\n");
        return 0;
    }

    len = strlen(path);
    _log_filename = calloc(len + 1, sizeof(char));

    if (!_log_filename) osh_malloc_err("osh_setup_logger() *_log_filename");
    memcpy(_log_filename, path, len);

    _fplog = fopen(_log_filename,"w");
    if (_fplog == NULL) {
        fprintf(stderr,"%scannot open >%s< for writing. %s\n", OSH_LOG_PREFIX_ERROR, _log_filename, strerror(errno));
        free(_log_filename);
        exit(EX_IOERR);
    }

    _log_level = log_level;
    _logger_save_active = 1;

    /* setup message database */
    _message_db_init();

    return 1;
}


int osh_close_logger(void) {

    _message_db_summary(); /* print any occuring warning messages */

    if (_fplog == NULL) {
        fprintf(stderr,"%scannot close %s for writing.\n", OSH_LOG_PREFIX_ERROR, _log_filename);
        exit(EX_IOERR);
    }
    fclose(_fplog);
    _logger_save_active = 0;
    free(_log_filename);
    return 0;
}


int osh_set_loglevel(int log_level){
    if (!_logger_save_active) {
        return 1;
    }
    _log_level = log_level;
    return 0;
}


int osh_get_loglevel(void){
    return _log_level;
}


void osh_err(int status, char const *msg, ...) {
    char s[_MAX_LINE_LENGTH];

    va_list args1, args2;
    va_start(args1, msg);
    va_copy(args2, args1);

    memset(s, 0, _MAX_LINE_LENGTH);
    snprintf(s, _MAX_LINE_LENGTH, "%s%s\n", OSH_LOG_PREFIX_ERROR, msg);

    if (_log_level >= OSH_LOG_ERR) {
        /* Print to stderr */
        vfprintf(stderr, s, args1);
    }

    /* all ERR will be saved to logger, irrespectly of log level, if save is active */
    if(_logger_save_active) {
        vfprintf(_fplog, s, args2);
        fflush(_fplog);
    }

    va_end(args1);
    va_end(args2);

    exit(status);
}


void osh_malloc_err(char const *msg, ...) {
    char s[_MAX_LINE_LENGTH];

    va_list args1, args2;
    va_start(args1, msg);
    va_copy(args2, args1);

    memset(s, 0, _MAX_LINE_LENGTH);
    snprintf(s, _MAX_LINE_LENGTH, "%s Could not allocate memory. malloc() %s\n", OSH_LOG_PREFIX_ERROR, msg);

    if (_log_level >= OSH_LOG_ERR) {
        vfprintf(stderr, s, args1);
    }

    /* all ERR will be saved to logger, irrespectly of log level, if save is active */
    if(_logger_save_active) {
        vfprintf(_fplog, s, args2);
        fflush(_fplog);
    }

    va_end(args1);
    va_end(args2);

    exit(EX_UNAVAILABLE);
}


void osh_warn(char const *msg, ...) {
    char s[_MAX_LINE_LENGTH];
    va_list args1, args2;
    va_start(args1, msg);
    va_copy(args2, args1);

    memset(s, 0, _MAX_LINE_LENGTH);
    snprintf(s, _MAX_LINE_LENGTH - 2, "%s%s", OSH_LOG_PREFIX_WARN, msg);

    /* check database */
    if(!_message_db_check(s, OSH_LOG_WARN)) {
        return;
    }

    if (_log_level >= OSH_LOG_WARN) {
        va_start(args1, msg);
        vfprintf(stderr, s, args1);
    }

    /* all warnings will be saved to logger, irrespectly of log level, if save is active */
    if(_logger_save_active) {
        va_start(args2, msg);
        vfprintf(_fplog, s, args2);
        fflush(_fplog);
    }

    va_end(args1);
    va_end(args2);
}


void osh_info(char const *msg, ...) {
    char s[_MAX_LINE_LENGTH] = "";

    va_list args1, args2;
    va_start(args1, msg);
    va_copy(args2, args1);

    memset(s, 0, _MAX_LINE_LENGTH);
    snprintf(s, _MAX_LINE_LENGTH - 2, "%s%s", OSH_LOG_PREFIX_INFO, msg);

    if (_log_level >= OSH_LOG_INFO) {
        va_start(args1, msg);
        vfprintf(stdout, s, args1);
        fflush(stdout);
    }

    /* all INFO will be saved to logger, irrespectly of log level, if save is active */
    if(_logger_save_active) {
        va_start(args2, msg);
        vfprintf(_fplog, s, args2);
        fflush(_fplog);
    }

    va_end(args1);
    va_end(args2);
}


void osh_debug(char const *msg, ...) {
    char s[_MAX_LINE_LENGTH];
    va_list args1, args2;
    va_start(args1, msg);
    va_copy(args2, args1);

    memset(s, 0, _MAX_LINE_LENGTH);
    snprintf(s, _MAX_LINE_LENGTH - 2, "%s%s", OSH_LOG_PREFIX_DEBUG, msg);

    if (_log_level >= OSH_LOG_DEBUG) {
        va_start(args1, msg);
        vfprintf(stderr, s, args1);

        /* DEBUG will be saved to logger, if log-level sufficient */
        if(_logger_save_active) {
            va_start(args2, msg);
            vfprintf(_fplog, s, args2);
            fflush(_fplog);
        }
    }

    va_end(args1);
    va_end(args2);
}


void osh_log(char const *msg, ...) {
    char s[_MAX_LINE_LENGTH];
    va_list args;

    if (_logger_save_active) {

        memset(s, 0, _MAX_LINE_LENGTH);
        snprintf(s, _MAX_LINE_LENGTH - 2, "%s", msg);

        va_start(args, msg);
        vfprintf(_fplog, s, args);
        fflush(_fplog);
        va_end(args);
    }
}


/**
 * @brief Add message to message database. If message has is known has been occuring more then OSH_LOG_REPEAT times, then return 0, else 1.
 * The database will also count how many time a message has been occuring.
 *
 * @param[in] s - message string
 * @param[in] level - message level (not really used for now, but may be used in future)
 *
 * @returns 1 if message is new, 0 if message has been occuring more then OSH_LOG_REPEAT times
 *
 * @author Niels Bassler
 */
int _message_db_check(char const *s, int level) {

    int message_index;
    int i;

    /* Search for the message in the database */
    for (i = 0; i < message_db.len; i++) {
        if (strcmp(message_db.msg[i], s) == 0) {
            /* Message found, increment count */
            message_db.count[i]++;
            /* Check if message has occurred more than OSH_LOG_REPEAT times */
            if (message_db.count[i] > OSH_LOG_REPEAT) {
                /* Message has occurred too many times, return 0 */
                return 0;
            } else {
                /* Message count is within limit, return 1 */
                if(message_db.count[i] == OSH_LOG_REPEAT) {
                    osh_info(
                        "The following warning message has been repeated %d times and further repetitions will be suppressed:\n",
                        OSH_LOG_REPEAT);
                }
                return 1;
            }
        }
    }

    /* Message not found, add new message to the database */
    message_index = message_db.len; /* new message index */
    message_db.len++;
    message_db.msg = realloc(message_db.msg, message_db.len * sizeof(char *));
    message_db.count = realloc(message_db.count, message_db.len * sizeof(int));
    message_db.level = realloc(message_db.level, message_db.len * sizeof(int));

    /*  Check for allocation errors */
    if (!message_db.msg || !message_db.count || !message_db.level) {
        osh_malloc_err("osh_message_db() extending database");
    }

    /* Allocate memory for the new message and copy it */
    message_db.msg[message_index] = strdup(s);
    if (!message_db.msg[message_index]) {
        osh_malloc_err("osh_message_db() message_db.msg[message_index]");
    }

    /* Initialize count and level for the new message */
    message_db.count[message_index] = 1;
    message_db.level[message_index] = level;

    /* New message added successfully, return 1 */
    return 1;
}


/**
 * @brief Print summary of the messages recorded in the message database, and the number of occurences.
 *
 * @author Niels Bassler
 */
void _message_db_summary(void) {
    int len;
    int i;

    if (message_db.len == 0) {
        return;
    }

    osh_info("\n");
    osh_info("Warning messages summary:\n");
    for (i = 0; i < message_db.len; i++) {

        /* Check if the last character is a newline */
        len = strlen(message_db.msg[i]);
        if (len > 0 && message_db.msg[i][len - 1] == '\n') {
            /* Temporarily remove the newline for printing */
            message_db.msg[i][len - 1] = '\0';
        }

        osh_info ("\"%s \" - count:%d \n",
                  message_db.msg[i],
                  message_db.count[i]);

    } /* end of for loop */
}


/* initialize message db */
void _message_db_init(void) {
    message_db.msg = NULL;
    message_db.count = NULL;
    message_db.len = 0;

    message_db.msg = calloc(1, sizeof(char *));
    if (!message_db.msg) osh_malloc_err("osh_message_db ()message_db.msg");
    message_db.count = calloc(1, sizeof(int));
    if (!message_db.count) osh_malloc_err("osh_message_db ()message_db.count");
    message_db.level = calloc(1, sizeof(int));
    if (!message_db.level) osh_malloc_err("osh_message_db ()message_db.level");

    message_db.count[0] = 0;

}


void _message_db_free(void) {
    int i;
    /* Free each individual message string */
    for (i = 0; i < message_db.len; i++) {
        free(message_db.msg[i]);
    }

    /* Free the arrays for messages, counts, and levels */
    free(message_db.msg);
    free(message_db.count);
    free(message_db.level);

    /* Reset the database structure fields */
    message_db.msg = NULL;
    message_db.count = NULL;
    message_db.level = NULL;
    message_db.len = 0;
}