#ifndef _OSH_LOGGER
#define _OSH_LOGGER

#include "osh_exit.h"

#define OSH_LOG_TRACE   6
#define OSH_LOG_DEBUG   5
#define OSH_LOG_INFO    4  /* default level */
#define OSH_LOG_WARN    3
#define OSH_LOG_ERR     2
#define OSH_LOG_FATAL   1
#define OSH_LOG_OFF     0

/* prefix will prepended any string used by the corresponding OSH_ log function */
#define OSH_LOG_PREFIX_ERROR   "*** Error: "
#define OSH_LOG_PREFIX_WARN    "*** Warning: "
#define OSH_LOG_PREFIX_INFO    "  "
#define OSH_LOG_PREFIX_DEBUG   "Debug: "
#define OSH_LOG_PREFIX_LOG     ""
#define OSH_LOG_HLINE          "-----------------------------------------------------\n"

#define OSH_LOG_REPEAT         10  /* how often a message is repeated before it is suppressed */


/**
 * @brief Setup the logger for OpenShieldHIT
 *
 * @param[in] path - path to logfile
 * @param[in] log_level - see defines in osh_logger.h
 *
 * @returns 1 if successful, 0 if logger was already setup
 *
 * @warning This is a temporary solution, since it is not thread-safe.
 *
 * @author Niels Bassler
 */
int osh_setup_logger(char const *path, int log_level);

/**
 * @brief Shutdown the logger for OpenShieldHIT.
 *
 * @returns 0 if successful
 *
 * @author Niels Bassler
 */
int osh_close_logger(void);

/**
 * @brief set the logger level
 *
 * @param log_level : logger level e.g. osh_LOG_WARN
 * @returns 0 if successful, 1 if unsuccessful
 *
 * @author Niels Bassler
 */
int osh_set_loglevel(int log_level);

/**
 * @brief get current logger level
 *
 * @returns current logging level
 *
 * @author Niels Bassler
 */
int osh_get_loglevel(void);


/**
 * @brief Print message to STDERR (and logfile if configured) and terminate.
 *
 * Appends a newline automatically. Accepts printf-style format arguments.
 *
 * @param[in] status Process exit code (use EX_* from osh_exit.h).
 * @param[in] msg    Null-terminated format string.
 * @param[in] ...    Optional printf-style arguments.
 *
 * @author Niels Bassler
 */
void osh_err(int status, char const *msg, ...);

/**
 * @brief Print message to STDERR and exit with EX_UNAVAILABLE.
 * Write also to the logfile, if it was setup.
 *
 * Message must be terminated with a null byte: //CHAR(0)
 * Message will automatically be terminated with a newline.
 *
 * @param[in] msg - Null-terminated format string.
 * @param[in] ... Optional format arguments
 *
 * @author Niels Bassler
 */
void osh_malloc_err(char const *msg,  ...);

/**
 * @brief Print a warning to stdout and logfile, if it was setup.
 *
 * Message must be terminated with a null byte: //CHAR(0)
 *
 * @param[in] msg - Null-terminated format string.
 * @param[in] ... Optional format arguments
 *
 * @author Niels Bassler
 */
void osh_warn(char const *msg, ...);

/**
 * @brief Write a message to stdout _and_ the logfile, if it was setup.
 *
 * Message must be terminated with a null byte: //CHAR(0)
 * @param[in] msg - Null-terminated format string.
 * @param[in] ... Optional format arguments
 *
 * @author Niels Bassler
 */
void osh_info(char const *msg, ...);

/**
 * @brief Write a message to stdout _and_ the logfile, if it was setup.
 *
 * Message must be terminated with a null byte: //CHAR(0)
 * @param[in] msg - Null-terminated format string.
 * @param[in] ... Optional format arguments
 *
 * @author Niels Bassler
 */
void osh_debug(char const *msg, ...);

/**
 * @brief Write a message to the logfile, if it was setup. Do not output to stdout.
 *
 * Message must be terminated with a null byte: //CHAR(0)
 * @param[in] msg - Null-terminated format string.
 * @param[in] ... Optional format arguments
 *
 * @author Niels Bassler
 */
void osh_log(char const *msg, ...);

#endif /* !_OSH_LOGGER */
