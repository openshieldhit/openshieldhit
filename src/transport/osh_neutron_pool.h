#ifndef OSH_NEUTRON_POOL_H
#define OSH_NEUTRON_POOL_H

#include <stddef.h>

/**
 * @file osh_neutron_pool.h
 * @brief Neutron secondary pool — infrastructure stub for future neutron transport.
 *
 * @details
 * Neutrons produced by nuclear reactions (abrasion, Fermi break-up) are not
 * transported through the CSDA charged-particle loop.  They are collected
 * here instead.  Currently only the count is stored; actual position/direction/
 * energy arrays will be added when neutron transport is implemented.
 */
struct osh_neutron_pool {
    size_t n_created; /**< Neutrons produced but not yet transported. */
};

#endif /* OSH_NEUTRON_POOL_H */
