#include "apps/osh/osh_app_osh.h"

#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_beam_parse.h"
#include "beam/osh_beam_spots.h"
#include "openshieldhit/file.h"
#include "openshieldhit/logger.h"
#include "openshieldhit/status.h"

int osh_beam_setup_from_path(char const *path, struct osh_logger *lg, struct osh_beam_workspace **wb_out) {
    int rc = OSH_OK;
    struct oshfile *sf = NULL;
    struct osh_beam_workspace *wb = NULL;
    char *spotlist_path = NULL;

    (void) lg;

    if (!path || !wb_out) {
        return OSH_EINVAL;
    }
    *wb_out = NULL;

    sf = osh_fopen(path);
    if (!sf) {
        return OSH_EIO;
    }

    rc = osh_beam_workspace_create(&wb);
    if (rc != OSH_OK) {
        osh_fclose(sf);
        return rc;
    }

    rc = osh_beam_parse(sf, wb, &spotlist_path);
    osh_fclose(sf);
    sf = NULL;
    if (rc != OSH_OK) {
        free(spotlist_path);
        osh_beam_workspace_free(wb);
        return rc;
    }

    if (spotlist_path) {
        osh_info("Loading spotlist before beam post-parse");
        rc = osh_beam_spotlist_load(wb, spotlist_path);
        if (rc != OSH_OK) {
            free(spotlist_path);
            osh_beam_workspace_free(wb);
            return rc;
        }
    }
    free(spotlist_path);

    rc = osh_beam_workspace_prepare(wb);
    if (rc != OSH_OK) {
        osh_beam_workspace_free(wb);
        return rc;
    }

    *wb_out = wb;
    return OSH_OK;
}
