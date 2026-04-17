#ifndef OSH_TRANSPORT_SCHEDULER_H
#define OSH_TRANSPORT_SCHEDULER_H

#include "openshieldhit/status.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Internal header: scheduler seam for future multi-pool transport orchestration.
 * Not part of the public API yet.
 */

enum osh_transport_family {
    OSH_TRANSPORT_FAMILY_ION = 0,
    OSH_TRANSPORT_FAMILY_NEUTRON = 1,
    OSH_TRANSPORT_FAMILY_PHOTON = 2,
    OSH_TRANSPORT_FAMILY_ELECTRON = 3,
    OSH_TRANSPORT_FAMILY_COUNT = 4
};

/*
 * Minimal scheduler state for the future outer transport loop.
 *
 * The eventual multi-pool driver will own one queue or pool per transport
 * family and update has_work[] as families generate secondaries for each
 * other.  For now this struct only tracks which families are enabled for a
 * run and which of them currently have work pending.
 */
struct osh_transport_scheduler {
    char enabled[OSH_TRANSPORT_FAMILY_COUNT];
    char has_work[OSH_TRANSPORT_FAMILY_COUNT];
};

/**
 * @brief Reset the transport scheduler to an empty state.
 *
 * @details
 * Clears all enabled and has_work flags.  The caller must then enable the
 * families that participate in the current run and mark the ones that hold
 * pending work.
 *
 * @param[out] scheduler  Scheduler state to reset.
 */
void osh_transport_scheduler_reset(struct osh_transport_scheduler *scheduler);

/**
 * @brief Enable one transport family in the scheduler.
 *
 * @param[in,out] scheduler  Scheduler to modify.
 * @param[in]     family     Transport family to enable.
 *
 * @returns OSH_OK on success, OSH_EINVAL if any argument is invalid.
 */
enum osh_status osh_transport_scheduler_enable(struct osh_transport_scheduler *scheduler,
                                               enum osh_transport_family family);

/**
 * @brief Mark whether one transport family currently has queued work.
 *
 * @details
 * This is a scheduler-only bookkeeping update.  The owning transport loop is
 * still responsible for allocating, filling, and draining the actual pools.
 *
 * @param[in,out] scheduler  Scheduler to modify.
 * @param[in]     family     Transport family to update.
 * @param[in]     has_work   1 if the family queue is non-empty, 0 if empty.
 */
void osh_transport_scheduler_set_has_work(struct osh_transport_scheduler *scheduler,
                                          enum osh_transport_family family,
                                          int has_work);

/**
 * @brief Select the next enabled family that has pending work.
 *
 * @details
 * The current stub uses a fixed priority order:
 *   ions -> neutrons -> photons -> electrons
 *
 * This is deterministic and sufficient for the current ion-only driver.
 * When multiple active pools exist, this function is the intended owner of
 * any later round-robin or cost-based scheduling policy.
 *
 * @param[in]  scheduler   Scheduler state.
 * @param[out] family_out  Receives the selected family on success.
 *
 * @returns 1 if a runnable family was found, 0 otherwise.
 */
int osh_transport_scheduler_next(struct osh_transport_scheduler const *scheduler,
                                 enum osh_transport_family *family_out);

/**
 * @brief Return a human-readable name for a transport family.
 *
 * @param[in] family  Family enum value.
 *
 * @returns Static name string, or "invalid" for an unknown enum value.
 */
char const *osh_transport_family_name(enum osh_transport_family family);

#ifdef __cplusplus
}
#endif

#endif /* OSH_TRANSPORT_SCHEDULER_H */
