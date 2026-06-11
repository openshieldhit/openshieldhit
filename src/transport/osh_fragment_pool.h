#ifndef OSH_FRAGMENT_POOL_H
#define OSH_FRAGMENT_POOL_H

#include <stddef.h>

/**
 * @file osh_fragment_pool.h
 * @brief Nuclear fragment pool — counter for unprocessed residual fragments.
 *
 * @details
 * Light prefragments (A <= 16) from abrasion are de-excited in-event by the
 * Fermi break-up stage.  What is counted here are the *unprocessed* residues:
 * heavy fragments outside the break-up domain (A > 16), nuclides with no open
 * decay channel that are not whitelisted transport products, and truncation
 * overflow.  These await a future evaporation/SMM stage; only the count is
 * stored for now.
 */
struct osh_fragment_pool {
    size_t n_created;      /**< Residual nuclear fragments created but not processed. */
    size_t n_sent_breakup; /**< Prefragments handed to the Fermi break-up stage.      */
    size_t n_breakup;      /**< Prefragments de-excited (≥ 1 product emitted).        */
};

#endif /* OSH_FRAGMENT_POOL_H */
