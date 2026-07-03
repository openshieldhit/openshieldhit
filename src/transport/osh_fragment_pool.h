#ifndef OSH_FRAGMENT_POOL_H
#define OSH_FRAGMENT_POOL_H

#include <stddef.h>

/**
 * @file osh_fragment_pool.h
 * @brief Nuclear fragment pool — diagnostic counters for heavy recoil fragments.
 *
 * @details
 * Light prefragments (A <= 16) from abrasion are de-excited in-event by the
 * Fermi break-up stage.  The heavy end-products it emits — the final stable
 * break-up residues (C, N, O, ... recoils) and any prefragments outside the
 * break-up domain (A > 16) — are handed to event_out->fragments[] and then
 * transported by the ion step as recoil ions, or point-deposited when below the
 * transport threshold (see osh_transport_ion_step.c, issue #179).  This struct
 * only holds diagnostic counts; the fragments themselves live in the event.
 */
struct osh_fragment_pool {
    size_t n_created;      /**< Heavy recoil fragments emitted (transported/point-deposited). */
    size_t n_sent_breakup; /**< Prefragments handed to the Fermi break-up stage.      */
    size_t n_breakup;      /**< Prefragments de-excited (≥ 1 product emitted).        */
};

#endif /* OSH_FRAGMENT_POOL_H */
