#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_app_osh.h"
#include "openshieldhit/scoring.h"
#include "openshieldhit/status.h"
#include "scoring/runtime/osh_scoring_compile.h"
#include "scoring/runtime/osh_scoring_postprocess.h"
#include "scoring/save/osh_scoring_save.h"
#include "scoring/save/osh_scoring_save_bdo2019_raw.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

#define DETECT_PATH "osh_scoring_save_detect.tmp"
#define ASCII_PATH "out_ascii.txt"
#define BDO_PATH "out_binary.bdo"
#define CYL_ASCII_PATH "out_cyl_ascii.txt"
#define CYL_BDO_PATH "out_cyl_binary.bdo"
#define BDO_TEST_MAX_PAYLOAD 4096u

static void read_file_bytes(char const *path, unsigned char *buf, size_t nbytes);
static int file_token_contains_text(char const *path, unsigned long long tag_id, char const *needle);
static size_t file_token_read_llints(char const *path, unsigned long long tag_id, long long int *values, size_t cap);
static size_t bdo_payload_size(struct osh_scoring_bdo2019_tag const *tag);
static void write_detect_file(char const *content);

static void test_save_bdo2019_with_dose_and_dlet(void) {
    char const *detect_text = "Geometry Mesh\n"
                              "    Name G\n"
                              "    X 0 1 1\n"
                              "    Y 0 1 1\n"
                              "    Z 0 2 2\n"
                              "\n"
                              "Output\n"
                              "    Filename out_dlet.bdo\n"
                              "    Geo G\n"
                              "    Quantity DOSE\n"
                              "    Quantity DLET\n"
                              "    Quantity TLET\n"
                              "    Quantity DQEFF\n"
                              "    Quantity TQEFF\n"
                              "    Quantity DAVGE\n"
                              "    Quantity TAVGE\n"
                              "    Quantity DBETA\n"
                              "    Quantity TBETA\n";
    struct osh_scoring_workspace *ws;
    struct osh_scoring_runtime rt;
    unsigned char bdo_head[6];
    enum osh_status rc;
    size_t i;

    write_detect_file(detect_text);

    ws = NULL;
    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_setup_from_path(DETECT_PATH, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(rt.npages == 9u);

    /* Seed non-zero values and wire data2 (two-pass pages) with matching weights */
    for (i = 0; i < rt.npages; ++i) {
        rt.pages[i].acc.data[0] = 2.0;
        rt.pages[i].acc.data[1] = 4.0;
        if (rt.pages[i].has_data2) {
            rt.pages[i].acc.data2[0] = 1.0;
            rt.pages[i].acc.data2[1] = 1.0;
        }
    }

    rc = osh_scoring_postprocess(&rt);
    ASSERT_TRUE(rc == OSH_OK);

    rc = osh_scoring_save(ws, &rt, 100u);
    ASSERT_TRUE(rc == OSH_OK);

    read_file_bytes("out_dlet.bdo", bdo_head, sizeof(bdo_head));
    ASSERT_TRUE(memcmp(bdo_head, "xSH12A", 6u) == 0);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(DETECT_PATH);
    remove("out_dlet.bdo");
}

static void test_save_bdo2019_dirtydose_uses_sh12a_page_types(void) {
    char const *detect_text = "Geometry Mesh\n"
                              "    Name G\n"
                              "    X 0 1 1\n"
                              "    Y 0 1 1\n"
                              "    Z 0 1 1\n"
                              "\n"
                              "Output\n"
                              "    Filename out_dirtydose.bdo\n"
                              "    FileFormat BDO2019\n"
                              "    Geo G\n"
                              "    Quantity DirtyDose\n"
                              "    Quantity DirtyDoseGy\n";
    struct osh_scoring_workspace *ws;
    struct osh_scoring_runtime rt;
    long long int page_types[2];
    enum osh_status rc;

    write_detect_file(detect_text);

    ws = NULL;
    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_setup_from_path(DETECT_PATH, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(rt.npages == 2u);

    rc = osh_scoring_save(ws, &rt, 1u);
    ASSERT_TRUE(rc == OSH_OK);

    ASSERT_TRUE(file_token_read_llints("out_dirtydose.bdo", OSHBDO_PAG_TYPE, page_types, 2u) == 2u);
    ASSERT_TRUE(page_types[0] == 64);
    ASSERT_TRUE(page_types[1] == 65);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(DETECT_PATH);
    remove("out_dirtydose.bdo");
}

static void test_save_cyl_ascii_and_bdo(void) {
    char const *detect_text = "Geometry Cyl\n"
                              "    Name G\n"
                              "    R 0 4 2\n"
                              "    Z 0 2 2\n"
                              "\n"
                              "Output\n"
                              "    Filename out_cyl_ascii.txt\n"
                              "    FileFormat ASCII\n"
                              "    Geo G\n"
                              "    Quantity ENERGY\n"
                              "\n"
                              "Output\n"
                              "    Filename out_cyl_binary.bdo\n"
                              "    FileFormat BDO2019\n"
                              "    Geo G\n"
                              "    Quantity ENERGY\n";
    struct osh_scoring_workspace *ws;
    struct osh_scoring_runtime rt;
    unsigned char bdo_head[8];
    FILE *fp;
    char line[512];
    int saw_cyl_header;
    int saw_axis_header;
    int saw_first_data_row;
    enum osh_status rc;

    write_detect_file(detect_text);

    ws = NULL;
    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_setup_from_path(DETECT_PATH, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(ws != NULL);
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(rt.noutputs == 2u);
    ASSERT_TRUE(rt.npages == 2u);

    rt.pages[0].acc.data[0] = 1.0;
    rt.pages[0].acc.data[1] = 2.0;
    rt.pages[0].acc.data[2] = 3.0;
    rt.pages[0].acc.data[3] = 4.0;
    rt.pages[1].acc.data[0] = 10.0;
    rt.pages[1].acc.data[1] = 20.0;
    rt.pages[1].acc.data[2] = 30.0;
    rt.pages[1].acc.data[3] = 40.0;

    rc = osh_scoring_save(ws, &rt, 5u);
    ASSERT_TRUE(rc == OSH_OK);

    fp = fopen(CYL_ASCII_PATH, "r");
    ASSERT_TRUE(fp != NULL);
    saw_cyl_header = 0;
    saw_axis_header = 0;
    saw_first_data_row = 0;
    while (fgets(line, sizeof(line), fp) != NULL) {
        if (strstr(line, "# DETECTOR OUTPUT CYL") != NULL) {
            saw_cyl_header = 1;
        }
        if (strstr(line, "# Z R ENERGY") != NULL) {
            saw_axis_header = 1;
        }
        if (strstr(line, "5.000000000000e-01 1.000000000000e+00 2.000000000000e-01") != NULL) {
            saw_first_data_row = 1;
        }
    }
    ASSERT_TRUE(saw_cyl_header);
    ASSERT_TRUE(saw_axis_header);
    ASSERT_TRUE(saw_first_data_row);
    ASSERT_TRUE(fclose(fp) == 0);

    read_file_bytes(CYL_BDO_PATH, bdo_head, sizeof(bdo_head));
    ASSERT_TRUE(memcmp(bdo_head, "xSH12A", 6u) == 0);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(DETECT_PATH);
    remove(CYL_ASCII_PATH);
    remove(CYL_BDO_PATH);
}

static void test_save_ascii_rejects_mixed_diff_layout(void) {
    char const *detect_text = "Geometry Mesh\n"
                              "    Name G\n"
                              "    X 0 1 1\n"
                              "    Y 0 1 1\n"
                              "    Z 0 1 1\n"
                              "\n"
                              "Output\n"
                              "    Filename out_bad_ascii.txt\n"
                              "    FileFormat ASCII\n"
                              "    Geo G\n"
                              "    Quantity FLUENCE\n"
                              "    Diff1 0 10 5\n"
                              "    Diff1Type LET\n"
                              "    Quantity FLUENCE\n"
                              "    Diff1 0 8 4\n"
                              "    Diff1Type QEFF\n";
    struct osh_scoring_workspace *ws;
    struct osh_scoring_runtime rt;
    enum osh_status rc;

    write_detect_file(detect_text);

    ws = NULL;
    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_setup_from_path(DETECT_PATH, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);

    /* ASCII has one table shape per Output, so it cannot safely combine pages
     * with different diff axes.  BDO carries page-local axis metadata and remains
     * the format for heterogeneous differential outputs. */
    rc = osh_scoring_save(ws, &rt, 1u);
    ASSERT_TRUE(rc == OSH_ENOTSUP);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(DETECT_PATH);
    remove("out_bad_ascii.txt");
}

static void test_save_bdo2019_diff_log_units(void) {
    char const *detect_text = "Geometry Mesh\n"
                              "    Name G\n"
                              "    X 0 1 1\n"
                              "    Y 0 1 1\n"
                              "    Z 0 1 1\n"
                              "\n"
                              "Output\n"
                              "    Filename out_diff_log.bdo\n"
                              "    FileFormat BDO2019\n"
                              "    Geo G\n"
                              "    Quantity Fluence\n"
                              "    Diff1 1 100 2 LOG\n"
                              "    Diff1Type EKIN\n";
    struct osh_scoring_workspace *ws;
    struct osh_scoring_runtime rt;
    enum osh_status rc;

    write_detect_file(detect_text);

    ws = NULL;
    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_setup_from_path(DETECT_PATH, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(rt.npages == 1u);
    ASSERT_TRUE(rt.pages[0].diff_log == 1);

    rt.pages[0].acc.data[0] = 9.0;
    rt.pages[0].acc.data[1] = 90.0;
    rc = osh_scoring_postprocess(&rt);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(rt.pages[0].acc.data[0] == 1.0);
    ASSERT_TRUE(rt.pages[0].acc.data[1] == 1.0);

    rc = osh_scoring_save(ws, &rt, 1u);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(file_token_contains_text("out_diff_log.bdo", OSHBDO_PAG_DATA_UNIT, "/cm^2/MeV"));
    ASSERT_TRUE(file_token_contains_text("out_diff_log.bdo", OSHBDO_PAG_DIF_UNITS, "/cm^2;MeV"));

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(DETECT_PATH);
    remove("out_diff_log.bdo");
}

static void test_save_bdo2019_diff_slash_unit_parentheses(void) {
    char const *detect_text = "Geometry Mesh\n"
                              "    Name G\n"
                              "    X 0 1 1\n"
                              "    Y 0 1 1\n"
                              "    Z 0 1 1\n"
                              "\n"
                              "Output\n"
                              "    Filename out_diff_let.bdo\n"
                              "    FileFormat BDO2019\n"
                              "    Geo G\n"
                              "    Quantity Fluence\n"
                              "    Diff1 0 10 2\n"
                              "    Diff1Type DEDX\n";
    struct osh_scoring_workspace *ws;
    struct osh_scoring_runtime rt;
    enum osh_status rc;

    write_detect_file(detect_text);

    ws = NULL;
    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_setup_from_path(DETECT_PATH, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(rt.npages == 1u);

    rt.pages[0].acc.data[0] = 5.0;
    rt.pages[0].acc.data[1] = 10.0;
    rc = osh_scoring_postprocess(&rt);
    ASSERT_TRUE(rc == OSH_OK);

    rc = osh_scoring_save(ws, &rt, 1u);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(file_token_contains_text("out_diff_let.bdo", OSHBDO_PAG_DATA_UNIT, "/cm^2/(MeV/cm)"));
    ASSERT_TRUE(file_token_contains_text("out_diff_let.bdo", OSHBDO_PAG_DIF_UNITS, "/cm^2;MeV/cm"));

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(DETECT_PATH);
    remove("out_diff_let.bdo");
}

/* Zone ASCII + BDO save paths.  Zone geometry needs app-level zone-name
 * resolution before compile; this test fills the resolved transport zone ids
 * directly (as the step/point unit tests do) to exercise the library save path
 * without a geo.dat workspace. */
static void test_save_zone_ascii_and_bdo(void) {
    char const *detect_text = "Geometry Zone\n"
                              "    Name Z\n"
                              "    Zone Entrance\n"
                              "    Volume 2.0\n"
                              "    Zone Target\n"
                              "    Volume 4.0\n"
                              "\n"
                              "Output\n"
                              "    Filename out_zone_ascii.txt\n"
                              "    FileFormat ASCII\n"
                              "    Geo Z\n"
                              "    Quantity Energy\n"
                              "    Quantity Fluence\n"
                              "    Quantity Dose\n"
                              "\n"
                              "Output\n"
                              "    Filename out_zone_binary.bdo\n"
                              "    FileFormat BDO2019\n"
                              "    Geo Z\n"
                              "    Quantity Energy\n"
                              "    Quantity Fluence\n"
                              "    Quantity Dose\n";
    struct osh_scoring_workspace *ws;
    struct osh_scoring_runtime rt;
    unsigned char bdo_head[8];
    FILE *fp;
    char line[512];
    int saw_zone_header;
    int saw_zone_bin;
    int saw_col_header;
    int saw_zone3_row;
    int saw_zone7_row;
    size_t i;
    enum osh_status rc;

    write_detect_file(detect_text);

    ws = NULL;
    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_setup_from_path(DETECT_PATH, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(ws != NULL);
    ASSERT_TRUE(ws->ngeometries == 1u);
    ASSERT_TRUE(ws->geometries[0].nzone_indices == 2u);

    /* Resolve the two Zone selectors to transport zone ids 3 and 7. */
    ws->geometries[0].zone_indices = (size_t *) calloc(2u, sizeof(*ws->geometries[0].zone_indices));
    ASSERT_TRUE(ws->geometries[0].zone_indices != NULL);
    ws->geometries[0].zone_indices[0] = 3u;
    ws->geometries[0].zone_indices[1] = 7u;

    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(rt.noutputs == 2u);
    ASSERT_TRUE(rt.npages == 6u);

    /* Seed the two zone bins of every page by score kind.  ENERGY is written
     * raw; FLUENCE/DOSE get ÷volume in postprocess (vol 2 and 4 cm3 -> ×0.5,
     * ×0.25). */
    for (i = 0; i < rt.npages; ++i) {
        switch (rt.pages[i].score_kind) {
        case OSH_SCORING_SCORE_ENERGY:
            rt.pages[i].acc.data[0] = 10.0;
            rt.pages[i].acc.data[1] = 20.0;
            break;
        case OSH_SCORING_SCORE_FLUENCE:
            rt.pages[i].acc.data[0] = 4.0; /* -> 4*0.5 = 2 */
            rt.pages[i].acc.data[1] = 8.0; /* -> 8*0.25 = 2 */
            break;
        case OSH_SCORING_SCORE_DOSE:
            rt.pages[i].acc.data[0] = 6.0;  /* -> 6*0.5 = 3 */
            rt.pages[i].acc.data[1] = 12.0; /* -> 12*0.25 = 3 */
            break;
        default:
            break;
        }
    }

    rc = osh_scoring_postprocess(&rt);
    ASSERT_TRUE(rc == OSH_OK);

    /* nstat = 10: NORM columns divide by nstat at write.  Per zone the rows are
     * ENERGY (1.0 / 2.0), FLUENCE (0.2 / 0.2), DOSE (0.3 / 0.3). */
    rc = osh_scoring_save(ws, &rt, 10u);
    ASSERT_TRUE(rc == OSH_OK);

    fp = fopen("out_zone_ascii.txt", "r");
    ASSERT_TRUE(fp != NULL);
    saw_zone_header = 0;
    saw_zone_bin = 0;
    saw_col_header = 0;
    saw_zone3_row = 0;
    saw_zone7_row = 0;
    while (fgets(line, sizeof(line), fp) != NULL) {
        if (strstr(line, "# DETECTOR OUTPUT ZONE") != NULL) {
            saw_zone_header = 1;
        }
        if (strstr(line, "# ZONE BIN:") != NULL) {
            saw_zone_bin = 1;
        }
        if (strstr(line, "# ZONE ENERGY FLUENCE DOSE") != NULL) {
            saw_col_header = 1;
        }
        /* First column is the resolved transport zone id, then the three NORM
         * columns divided by nstat. */
        if (strstr(line, "3 1.000000000000e+00 2.000000000000e-01 3.000000000000e-01") != NULL) {
            saw_zone3_row = 1;
        }
        if (strstr(line, "7 2.000000000000e+00 2.000000000000e-01 3.000000000000e-01") != NULL) {
            saw_zone7_row = 1;
        }
    }
    ASSERT_TRUE(saw_zone_header);
    ASSERT_TRUE(saw_zone_bin);
    ASSERT_TRUE(saw_col_header);
    ASSERT_TRUE(saw_zone3_row);
    ASSERT_TRUE(saw_zone7_row);
    ASSERT_TRUE(fclose(fp) == 0);

    read_file_bytes("out_zone_binary.bdo", bdo_head, sizeof(bdo_head));
    ASSERT_TRUE(memcmp(bdo_head, "xSH12A", 6u) == 0);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(DETECT_PATH);
    remove("out_zone_ascii.txt");
    remove("out_zone_binary.bdo");
}

int main(void) {
    char const *detect_text = "Geometry Mesh\n"
                              "    Name G\n"
                              "    X 0 2 2\n"
                              "    Y 0 1 1\n"
                              "    Z 0 2 2\n"
                              "\n"
                              "Output\n"
                              "    Filename out_ascii.txt\n"
                              "    FileFormat ASCII\n"
                              "    Geo G\n"
                              "    Quantity ENERGY\n"
                              "\n"
                              "Output\n"
                              "    Filename out_binary.bdo\n"
                              "    FileFormat BDO2019\n"
                              "    Geo G\n"
                              "    Quantity ENERGY\n";
    struct osh_scoring_workspace *ws;
    struct osh_scoring_runtime rt;
    unsigned char bdo_head[8];
    FILE *fp;
    char line[512];
    enum osh_status rc;

    write_detect_file(detect_text);

    ws = NULL;
    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_setup_from_path(DETECT_PATH, NULL, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(ws != NULL);
    rc = osh_scoring_compile(ws, NULL, &rt);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(rt.noutputs == 2u);
    ASSERT_TRUE(rt.npages == 2u);

    rt.pages[0].acc.data[0] = 1.0;
    rt.pages[0].acc.data[1] = 2.0;
    rt.pages[0].acc.data[2] = 3.0;
    rt.pages[0].acc.data[3] = 4.0;
    rt.pages[1].acc.data[0] = 10.0;
    rt.pages[1].acc.data[1] = 20.0;
    rt.pages[1].acc.data[2] = 30.0;
    rt.pages[1].acc.data[3] = 40.0;

    rc = osh_scoring_save(ws, &rt, 5u);
    ASSERT_TRUE(rc == OSH_OK);

    fp = fopen(ASCII_PATH, "r");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fgets(line, sizeof(line), fp) != NULL);
    ASSERT_TRUE(strstr(line, "OpenShieldHIT version") != NULL);
    while (fgets(line, sizeof(line), fp) != NULL) {
        /* data[0]=1.0 / nstat=5 = 0.2 */
        if (strstr(line, "5.000000000000e-01 5.000000000000e-01 5.000000000000e-01 2.000000000000e-01") != NULL) {
            break;
        }
    }
    ASSERT_TRUE(!feof(fp));
    ASSERT_TRUE(fclose(fp) == 0);

    read_file_bytes(BDO_PATH, bdo_head, sizeof(bdo_head));
    ASSERT_TRUE(memcmp(bdo_head, "xSH12A", 6u) == 0);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(DETECT_PATH);
    remove(ASCII_PATH);
    remove(BDO_PATH);

    test_save_bdo2019_with_dose_and_dlet();
    test_save_bdo2019_dirtydose_uses_sh12a_page_types();
    test_save_cyl_ascii_and_bdo();
    test_save_ascii_rejects_mixed_diff_layout();
    test_save_bdo2019_diff_log_units();
    test_save_bdo2019_diff_slash_unit_parentheses();
    test_save_zone_ascii_and_bdo();
    return 0;
}

static void read_file_bytes(char const *path, unsigned char *buf, size_t nbytes) {
    FILE *fp;

    fp = fopen(path, "rb");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fread(buf, 1u, nbytes, fp) == nbytes);
    ASSERT_TRUE(fclose(fp) == 0);
}

static int file_token_contains_text(char const *path, unsigned long long tag_id, char const *needle) {
    FILE *fp;
    struct osh_scoring_bdo2019_tag tag;
    size_t payload_size;
    char payload[BDO_TEST_MAX_PAYLOAD];
    size_t needle_len;
    size_t i;
    int found;

    fp = fopen(path, "rb");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fseek(fp, 24L, SEEK_SET) == 0); /* magic + endian + version field */
    needle_len = strlen(needle);
    found = 0;
    while (fread(&tag, sizeof(tag), 1u, fp) == 1u) {
        payload_size = bdo_payload_size(&tag);
        if (payload_size > sizeof(payload)) {
            break;
        }
        if (payload_size > 0u) {
            ASSERT_TRUE(fread(payload, 1u, payload_size, fp) == payload_size);
        }
        if (tag.tag == tag_id && needle_len <= payload_size) {
            for (i = 0u; i <= payload_size - needle_len; ++i) {
                if (memcmp(payload + i, needle, needle_len) == 0) {
                    found = 1;
                    break;
                }
            }
        }
        if (found) {
            break;
        }
    }
    ASSERT_TRUE(fclose(fp) == 0);
    return found;
}

static size_t file_token_read_llints(char const *path, unsigned long long tag_id, long long int *values, size_t cap) {
    FILE *fp;
    struct osh_scoring_bdo2019_tag tag;
    size_t payload_size;
    size_t nread;
    unsigned char payload[BDO_TEST_MAX_PAYLOAD];

    fp = fopen(path, "rb");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fseek(fp, 24L, SEEK_SET) == 0); /* magic + endian + version field */
    nread = 0u;
    while (fread(&tag, sizeof(tag), 1u, fp) == 1u) {
        payload_size = bdo_payload_size(&tag);
        ASSERT_TRUE(payload_size <= sizeof(payload));
        if (payload_size > 0u) {
            ASSERT_TRUE(fread(payload, 1u, payload_size, fp) == payload_size);
        }
        if (tag.tag == tag_id && strstr(tag.pltype, "i8") && payload_size == 8u && nread < cap) {
            memcpy(&values[nread], payload, sizeof(values[nread]));
            ++nread;
        }
    }
    ASSERT_TRUE(fclose(fp) == 0);
    return nread;
}

static size_t bdo_payload_size(struct osh_scoring_bdo2019_tag const *tag) {
    if (!tag) {
        return 0u;
    }
    if (tag->pltype[0] == 'S') {
        return (size_t) strtoull(tag->pltype + 1u, NULL, 10) * (size_t) tag->len;
    }
    if (strstr(tag->pltype, "f8") || strstr(tag->pltype, "i8") || strstr(tag->pltype, "u8")) {
        return 8u * (size_t) tag->len;
    }
    return 0u;
}

static void write_detect_file(char const *content) {
    FILE *fp;

    fp = fopen(DETECT_PATH, "w");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fputs(content, fp) >= 0);
    ASSERT_TRUE(fclose(fp) == 0);
}
