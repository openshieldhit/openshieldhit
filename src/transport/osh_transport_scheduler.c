#include "transport/osh_transport_scheduler.h"

#include <stddef.h>

/**
 * @brief Reset the transport scheduler to an empty state.
 *
 * @param[out] scheduler  Scheduler state to reset.
 */
void osh_transport_scheduler_reset(struct osh_transport_scheduler *scheduler) {
    size_t i;

    if (!scheduler) {
        return;
    }

    for (i = 0u; i < (size_t) OSH_TRANSPORT_FAMILY_COUNT; ++i) {
        scheduler->enabled[i] = 0;
        scheduler->has_work[i] = 0;
    }
}

/**
 * @brief Enable one transport family in the scheduler.
 *
 * @param[in,out] scheduler  Scheduler to modify.
 * @param[in]     family     Transport family to enable.
 *
 * @returns OSH_OK on success, OSH_EINVAL if any argument is invalid.
 */
enum osh_status osh_transport_scheduler_enable(struct osh_transport_scheduler *scheduler,
                                               enum osh_transport_family family) {
    if (!scheduler || family < OSH_TRANSPORT_FAMILY_ION || family >= OSH_TRANSPORT_FAMILY_COUNT) {
        return OSH_EINVAL;
    }

    scheduler->enabled[(size_t) family] = 1;
    return OSH_OK;
}

/**
 * @brief Mark whether one transport family currently has queued work.
 *
 * @param[in,out] scheduler  Scheduler to modify.
 * @param[in]     family     Transport family to update.
 * @param[in]     has_work   1 if the family queue is non-empty, 0 if empty.
 */
void osh_transport_scheduler_set_has_work(struct osh_transport_scheduler *scheduler,
                                          enum osh_transport_family family,
                                          int has_work) {
    if (!scheduler || family < OSH_TRANSPORT_FAMILY_ION || family >= OSH_TRANSPORT_FAMILY_COUNT) {
        return;
    }

    scheduler->has_work[(size_t) family] = has_work ? 1 : 0;
}

/**
 * @brief Select the next enabled family that has pending work.
 *
 * @param[in]  scheduler   Scheduler state.
 * @param[out] family_out  Receives the selected family on success.
 *
 * @returns 1 if a runnable family was found, 0 otherwise.
 */
int osh_transport_scheduler_next(struct osh_transport_scheduler const *scheduler,
                                 enum osh_transport_family *family_out) {
    size_t i;

    if (!scheduler || !family_out) {
        return 0;
    }

    for (i = 0u; i < (size_t) OSH_TRANSPORT_FAMILY_COUNT; ++i) {
        if (scheduler->enabled[i] && scheduler->has_work[i]) {
            *family_out = (enum osh_transport_family) i;
            return 1;
        }
    }

    return 0;
}

/**
 * @brief Return a human-readable name for a transport family.
 *
 * @param[in] family  Family enum value.
 *
 * @returns Static name string, or "invalid" for an unknown enum value.
 */
char const *osh_transport_family_name(enum osh_transport_family family) {
    switch (family) {
    case OSH_TRANSPORT_FAMILY_ION:
        return "ion";
    case OSH_TRANSPORT_FAMILY_NEUTRON:
        return "neutron";
    case OSH_TRANSPORT_FAMILY_PHOTON:
        return "photon";
    case OSH_TRANSPORT_FAMILY_ELECTRON:
        return "electron";
    case OSH_TRANSPORT_FAMILY_COUNT:
        break;
    }

    return "invalid";
}
