#include "scoring/save/osh_scoring_save_bdo2019_raw.h"

#include <stdlib.h>
#include <string.h>

static enum osh_status write_bytes(FILE *fp, void const *data, size_t size, size_t count);
static int host_is_big_endian(void);
static char const *host_endian_prefix(void);
static char const *host_dtype_prefix(void);
static void init_tag(struct osh_scoring_bdo2019_tag *tag, uint64_t tag_id, char const *pltype, uint64_t len);
static size_t padded_string_len(char const *str);

enum osh_status osh_scoring_bdo2019_write_preamble(FILE *fp, char const *version_string) {
    char version_field[16];
    enum osh_status rc;

    if (!fp || !version_string) {
        return OSH_EINVAL;
    }

    rc = write_bytes(fp, OSH_SCORING_BDO2019_MAGIC_NUMBER, 1u, 6u);
    if (rc != OSH_OK) {
        return rc;
    }
    rc = write_bytes(fp, host_endian_prefix(), 1u, 2u);
    if (rc != OSH_OK) {
        return rc;
    }

    memset(version_field, 0, sizeof(version_field));
    snprintf(version_field, sizeof(version_field), "%s", version_string);
    return write_bytes(fp, version_field, 1u, sizeof(version_field));
}

enum osh_status osh_scoring_bdo2019_write_token_str(FILE *fp, uint64_t tag_id, char const *str) {
    struct osh_scoring_bdo2019_tag tag;
    char *payload;
    size_t len;
    int nchar;
    enum osh_status rc;

    if (!fp || !str) {
        return OSH_EINVAL;
    }

    len = padded_string_len(str);
    payload = (char *) calloc(len, sizeof(*payload));
    if (!payload) {
        return OSH_ENOMEM;
    }
    memcpy(payload, str, strlen(str));

    init_tag(&tag, tag_id, "", 1u);
    nchar = snprintf(tag.pltype, sizeof(tag.pltype), "%s%zu", OSH_SCORING_BDO2019_PL_TYPE_CHAR, len);
    if (nchar < 0 || (size_t) nchar >= sizeof(tag.pltype)) {
        free(payload);
        return OSH_ESTATE;
    }

    rc = write_bytes(fp, &tag, sizeof(tag), 1u);
    if (rc == OSH_OK) {
        rc = write_bytes(fp, payload, 1u, len);
    }

    free(payload);
    return rc;
}

enum osh_status
osh_scoring_bdo2019_write_token_llint(FILE *fp, uint64_t tag_id, long long int const *values, size_t nvalues) {
    struct osh_scoring_bdo2019_tag tag;
    enum osh_status rc;
    int nchar;

    if (!fp || (!values && nvalues > 0u)) {
        return OSH_EINVAL;
    }

    init_tag(&tag, tag_id, "", (uint64_t) nvalues);
    nchar = snprintf(tag.pltype, sizeof(tag.pltype), "%s%s", host_dtype_prefix(), OSH_SCORING_BDO2019_PL_TYPE_LLSINT);
    if (nchar < 0 || (size_t) nchar >= sizeof(tag.pltype)) {
        return OSH_ESTATE;
    }

    rc = write_bytes(fp, &tag, sizeof(tag), 1u);
    if (rc != OSH_OK || nvalues == 0u) {
        return rc;
    }
    return write_bytes(fp, values, sizeof(*values), nvalues);
}

enum osh_status osh_scoring_bdo2019_write_token_int(FILE *fp, uint64_t tag_id, int const *values, size_t nvalues) {
    enum osh_status rc;
    long long int *tmp;
    size_t i;

    if (!fp || (!values && nvalues > 0u)) {
        return OSH_EINVAL;
    }

    if (nvalues == 0u) {
        return osh_scoring_bdo2019_write_token_llint(fp, tag_id, NULL, 0u);
    }

    tmp = (long long int *) calloc(nvalues, sizeof(*tmp));
    if (!tmp) {
        return OSH_ENOMEM;
    }
    for (i = 0; i < nvalues; ++i) {
        tmp[i] = (long long int) values[i];
    }

    rc = osh_scoring_bdo2019_write_token_llint(fp, tag_id, tmp, nvalues);
    free(tmp);
    return rc;
}

enum osh_status
osh_scoring_bdo2019_write_token_double(FILE *fp, uint64_t tag_id, double const *values, size_t nvalues) {
    struct osh_scoring_bdo2019_tag tag;
    enum osh_status rc;
    int nchar;

    if (!fp || (!values && nvalues > 0u)) {
        return OSH_EINVAL;
    }

    init_tag(&tag, tag_id, "", (uint64_t) nvalues);
    nchar = snprintf(tag.pltype, sizeof(tag.pltype), "%s%s", host_dtype_prefix(), OSH_SCORING_BDO2019_PL_TYPE_DOUBLE);
    if (nchar < 0 || (size_t) nchar >= sizeof(tag.pltype)) {
        return OSH_ESTATE;
    }

    rc = write_bytes(fp, &tag, sizeof(tag), 1u);
    if (rc != OSH_OK || nvalues == 0u) {
        return rc;
    }
    return write_bytes(fp, values, sizeof(*values), nvalues);
}

enum osh_status osh_scoring_bdo2019_write_token_float(FILE *fp, uint64_t tag_id, float const *values, size_t nvalues) {
    enum osh_status rc;
    double *tmp;
    size_t i;

    if (!fp || (!values && nvalues > 0u)) {
        return OSH_EINVAL;
    }

    if (nvalues == 0u) {
        return osh_scoring_bdo2019_write_token_double(fp, tag_id, NULL, 0u);
    }

    tmp = (double *) calloc(nvalues, sizeof(*tmp));
    if (!tmp) {
        return OSH_ENOMEM;
    }
    for (i = 0; i < nvalues; ++i) {
        tmp[i] = (double) values[i];
    }

    rc = osh_scoring_bdo2019_write_token_double(fp, tag_id, tmp, nvalues);
    free(tmp);
    return rc;
}

static enum osh_status write_bytes(FILE *fp, void const *data, size_t size, size_t count) {
    if (count == 0u) {
        return OSH_OK;
    }
    if (fwrite(data, size, count, fp) != count) {
        return OSH_EIO;
    }
    return OSH_OK;
}

static int host_is_big_endian(void) {
    uint16_t value = 0x0102u;
    unsigned char const *bytes = (unsigned char const *) &value;
    return bytes[0] == 0x01u;
}

static char const *host_endian_prefix(void) {
    return host_is_big_endian() ? OSH_SCORING_BDO2019_ENDIAN_BIG : OSH_SCORING_BDO2019_ENDIAN_LITTLE;
}

static char const *host_dtype_prefix(void) {
    return host_is_big_endian() ? OSH_SCORING_BDO2019_DTYPE_ENDIAN_BIG : OSH_SCORING_BDO2019_DTYPE_ENDIAN_LITTLE;
}

static void init_tag(struct osh_scoring_bdo2019_tag *tag, uint64_t tag_id, char const *pltype, uint64_t len) {
    memset(tag, 0, sizeof(*tag));
    tag->tag = tag_id;
    tag->len = len;
    if (pltype && pltype[0] != '\0') {
        snprintf(tag->pltype, sizeof(tag->pltype), "%s", pltype);
    }
}

static size_t padded_string_len(char const *str) {
    size_t len = strlen(str) + 1u;
    size_t rem = len % 8u;

    if (rem == 0u) {
        return len;
    }
    return len + (8u - rem);
}
