#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_app_osh.h"
#include "openshieldhit/scoring.h"
#include "openshieldhit/status.h"
#include "scoring/runtime/osh_scoring_compile.h"
#include "scoring/runtime/osh_scoring_postprocess.h"
#include "scoring/save/osh_scoring_save.h"

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

static void read_file_bytes(char const *path, unsigned char *buf, size_t nbytes);
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
                              "    Quantity TQEFF\n";
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
    ASSERT_TRUE(rt.npages == 5u);

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
    test_save_cyl_ascii_and_bdo();
    return 0;
}

static void read_file_bytes(char const *path, unsigned char *buf, size_t nbytes) {
    FILE *fp;

    fp = fopen(path, "rb");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fread(buf, 1u, nbytes, fp) == nbytes);
    ASSERT_TRUE(fclose(fp) == 0);
}

static void write_detect_file(char const *content) {
    FILE *fp;

    fp = fopen(DETECT_PATH, "w");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fputs(content, fp) >= 0);
    ASSERT_TRUE(fclose(fp) == 0);
}
