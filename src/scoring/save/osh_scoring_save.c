#include "scoring/save/osh_scoring_save.h"

#include <string.h>

#include "scoring/save/osh_scoring_save_ascii.h"
#include "scoring/save/osh_scoring_save_bdo2019.h"

static enum osh_status save_one_output(struct osh_scoring_save_request const *req, size_t output_idx);
static int fileformat_is_ascii(char const *fileformat);
static int fileformat_is_bdo2019(char const *fileformat);

enum osh_status osh_scoring_save(struct osh_scoring_save_request const *req) {
    enum osh_status rc;
    size_t i;

    if (!req || !req->ws || !req->rt) {
        return OSH_EINVAL;
    }
    if (req->ws->noutputs != req->rt->noutputs) {
        return OSH_ESTATE;
    }

    /* TODO: once asynchronous saving is added, this loop is the natural
     * hand-off point to a CPU-side writer queue. Keep runtime ownership
     * explicit and avoid mixing this path with transport execution state. */
    for (i = 0; i < req->rt->noutputs; ++i) {
        rc = save_one_output(req, i);
        if (rc != OSH_OK) {
            return rc;
        }
    }

    return OSH_OK;
}

static enum osh_status save_one_output(struct osh_scoring_save_request const *req, size_t output_idx) {
    char const *fileformat;

    if (output_idx >= req->ws->noutputs || output_idx >= req->rt->noutputs) {
        return OSH_EINVAL;
    }

    fileformat = req->ws->outputs[output_idx].fileformat;
    if (fileformat_is_ascii(fileformat)) {
        return osh_scoring_save_ascii_output(req, output_idx);
    }
    if (fileformat_is_bdo2019(fileformat)) {
        return osh_scoring_save_bdo2019_output(req, output_idx);
    }

    return OSH_ENOTSUP;
}

static int fileformat_is_ascii(char const *fileformat) {
    if (!fileformat) {
        return 0;
    }
    return (strcmp(fileformat, "text") == 0) || (strcmp(fileformat, "txt") == 0) || (strcmp(fileformat, "ascii") == 0)
           || (strcmp(fileformat, "dat") == 0);
}

static int fileformat_is_bdo2019(char const *fileformat) {
    if (!fileformat) {
        return 1;
    }
    return (strcmp(fileformat, "bdo") == 0) || (strcmp(fileformat, "bdo2019") == 0)
           || (strcmp(fileformat, "binary") == 0) || (strcmp(fileformat, "bin") == 0);
}
