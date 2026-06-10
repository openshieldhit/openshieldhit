#ifndef OSH_FRAGMENT_POOL_H
#define OSH_FRAGMENT_POOL_H

#include <stddef.h>

/**
 * @file osh_fragment_pool.h
 * @brief Nuclear fragment pool — infrastructure stub for future breakup transport.
 *
 * @details
 * Residual nuclear fragments produced by abrasion are collected here for the
 * future Fermi-breakup/SMM stage.  Currently only the count is stored; fragment
 * identity, excitation, position, and decay products will be added when that
 * physics is wired in.
 */
struct osh_fragment_pool {
    size_t n_created; /**< Residual nuclear fragments created but not processed. */
};

#endif /* OSH_FRAGMENT_POOL_H */
