#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "scoring/save/osh_scoring_save_bdo2019_raw.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

static FILE *osh_test_tmpfile(void);
static char const *host_endian_prefix(void);

static void test_tag_layout_is_24_bytes(void) {
    ASSERT_TRUE(sizeof(struct osh_scoring_bdo2019_tag) == 24u);
}

static void test_write_preamble_writes_magic_endian_and_version(void) {
    FILE *fp;
    char buf[24];
    size_t nread;

    fp = osh_test_tmpfile();
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(osh_scoring_bdo2019_write_preamble(fp, "v1.2.3") == OSH_OK);

    ASSERT_TRUE(fseek(fp, 0, SEEK_SET) == 0);
    nread = fread(buf, 1u, sizeof(buf), fp);
    ASSERT_TRUE(nread == sizeof(buf));

    ASSERT_TRUE(memcmp(buf, OSH_SCORING_BDO2019_MAGIC_NUMBER, 6u) == 0);
    ASSERT_TRUE(memcmp(buf + 6, host_endian_prefix(), 2u) == 0);
    ASSERT_TRUE(strcmp(buf + 8, "v1.2.3") == 0);

    ASSERT_TRUE(fclose(fp) == 0);
}

static void test_write_token_str_pads_payload_to_eight_bytes(void) {
    FILE *fp;
    struct osh_scoring_bdo2019_tag tag;
    char payload[8];

    fp = osh_test_tmpfile();
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(osh_scoring_bdo2019_write_token_str(fp, OSHBDO_SHVERSION, "abc") == OSH_OK);

    ASSERT_TRUE(fseek(fp, 0, SEEK_SET) == 0);
    ASSERT_TRUE(fread(&tag, sizeof(tag), 1u, fp) == 1u);
    ASSERT_TRUE(fread(payload, sizeof(payload), 1u, fp) == 1u);

    ASSERT_TRUE(tag.tag == (uint64_t) OSHBDO_SHVERSION);
    ASSERT_TRUE(tag.len == 1u);
    ASSERT_TRUE(strcmp(tag.pltype, "S8") == 0);
    ASSERT_TRUE(memcmp(payload, "abc\0\0\0\0\0", 8u) == 0);

    ASSERT_TRUE(fclose(fp) == 0);
}

static void test_write_token_int_writes_all_values_as_signed_64bit(void) {
    FILE *fp;
    struct osh_scoring_bdo2019_tag tag;
    int values[3] = {1, -2, 7};
    long long int payload[3];
    char pltype[8];

    fp = osh_test_tmpfile();
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(osh_scoring_bdo2019_write_token_int(fp, OSHBDO_EST_NPAGES, values, 3u) == OSH_OK);

    ASSERT_TRUE(fseek(fp, 0, SEEK_SET) == 0);
    ASSERT_TRUE(fread(&tag, sizeof(tag), 1u, fp) == 1u);
    ASSERT_TRUE(fread(payload, sizeof(payload), 1u, fp) == 1u);

    snprintf(pltype, sizeof(pltype), "%s%s", host_endian_prefix(), OSH_SCORING_BDO2019_PL_TYPE_LLSINT);
    ASSERT_TRUE(tag.tag == (uint64_t) OSHBDO_EST_NPAGES);
    ASSERT_TRUE(tag.len == 3u);
    ASSERT_TRUE(strcmp(tag.pltype, pltype) == 0);
    ASSERT_TRUE(payload[0] == 1);
    ASSERT_TRUE(payload[1] == -2);
    ASSERT_TRUE(payload[2] == 7);

    ASSERT_TRUE(fclose(fp) == 0);
}

static FILE *osh_test_tmpfile(void) {
#if defined(_MSC_VER)
    FILE *fp = NULL;
    if (tmpfile_s(&fp) != 0) {
        return NULL;
    }
    return fp;
#else
    return tmpfile();
#endif
}

static char const *host_endian_prefix(void) {
    uint16_t value = 0x0102u;
    unsigned char const *bytes = (unsigned char const *) &value;
    return (bytes[0] == 0x01u) ? OSH_SCORING_BDO2019_ENDIAN_BIG : OSH_SCORING_BDO2019_ENDIAN_LITTLE;
}

int main(void) {
    test_tag_layout_is_24_bytes();
    test_write_preamble_writes_magic_endian_and_version();
    test_write_token_str_pads_payload_to_eight_bytes();
    test_write_token_int_writes_all_values_as_signed_64bit();
    return 0;
}
