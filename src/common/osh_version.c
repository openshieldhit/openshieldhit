#include "common/osh_version.h"

#include "openshieldhit/version.h"

char const *osh_version_string(void) {
    return OSH_VERSION;
}

int osh_version_major(void) {
    return OSH_VERSION_MAJOR;
}

int osh_version_minor(void) {
    return OSH_VERSION_MINOR;
}

int osh_version_patch(void) {
    return OSH_VERSION_PATCH;
}
