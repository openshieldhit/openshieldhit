#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "common/osh_rc.h"
#include "scoring/osh_scoring.h"
#include "scoring/runtime/osh_scoring_prepare.h"
#include "scoring/save/osh_scoring_save.h"

#define ASSERT_TRUE(cond)                                                                                              \
    do {                                                                                                               \
        if (!(cond)) {                                                                                                 \
            fprintf(stderr, "ASSERT FAILED: %s (%s:%d)\n", #cond, __FILE__, __LINE__);                                 \
            exit(1);                                                                                                   \
        }                                                                                                              \
    } while (0)

static void build_temp_dir(char *path);
static void read_file_bytes(char const *path, unsigned char *buf, size_t nbytes);
static void write_temp_file(char *path, size_t cap, char const *content);

int main(void) {
    char detect_path[512];
    char out_dir[] = "/tmp/osh_scoring_save_testXXXXXX";
    char ascii_path[512];
    char bdo_path[512];
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
    struct osh_scoring_save_request req;
    unsigned char bdo_head[8];
    FILE *fp;
    char line[512];
    enum osh_status rc;

    build_temp_dir(out_dir);
    write_temp_file(detect_path, sizeof(detect_path), detect_text);

    ws = NULL;
    memset(&rt, 0, sizeof(rt));
    rc = osh_scoring_setup_from_path(detect_path, &ws);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(ws != NULL);
    rc = osh_scoring_prepare(ws, &rt);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(rt.noutputs == 2u);
    ASSERT_TRUE(rt.npages == 2u);

    rt.pages[0].data[0] = 1.0;
    rt.pages[0].data[1] = 2.0;
    rt.pages[0].data[2] = 3.0;
    rt.pages[0].data[3] = 4.0;
    rt.pages[1].data[0] = 10.0;
    rt.pages[1].data[1] = 20.0;
    rt.pages[1].data[2] = 30.0;
    rt.pages[1].data[3] = 40.0;

    memset(&req, 0, sizeof(req));
    req.out_dir = out_dir;
    req.ws = ws;
    req.rt = &rt;
    req.nstat = 5u;
    req.has_nstat = 1;

    rc = osh_scoring_save(&req);
    ASSERT_TRUE(rc == OSH_OK);

    snprintf(ascii_path, sizeof(ascii_path), "%s/out_ascii.txt", out_dir);
    snprintf(bdo_path, sizeof(bdo_path), "%s/out_binary.bdo", out_dir);

    fp = fopen(ascii_path, "r");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fgets(line, sizeof(line), fp) != NULL);
    ASSERT_TRUE(strstr(line, "OpenShieldHIT version") != NULL);
    while (fgets(line, sizeof(line), fp) != NULL) {
        if (strstr(line, "5.000000000000e-01 5.000000000000e-01 5.000000000000e-01 1.000000000000e+00") != NULL) {
            break;
        }
    }
    ASSERT_TRUE(!feof(fp));
    ASSERT_TRUE(fclose(fp) == 0);

    read_file_bytes(bdo_path, bdo_head, sizeof(bdo_head));
    ASSERT_TRUE(memcmp(bdo_head, "xSH12A", 6u) == 0);

    osh_scoring_runtime_free(&rt);
    osh_scoring_workspace_free(ws);
    remove(detect_path);
    remove(ascii_path);
    remove(bdo_path);
    rmdir(out_dir);
    return 0;
}

static void build_temp_dir(char *path) {
    ASSERT_TRUE(mkdtemp(path) != NULL);
}

static void read_file_bytes(char const *path, unsigned char *buf, size_t nbytes) {
    FILE *fp;

    fp = fopen(path, "rb");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fread(buf, 1u, nbytes, fp) == nbytes);
    ASSERT_TRUE(fclose(fp) == 0);
}

static void write_temp_file(char *path, size_t cap, char const *content) {
    FILE *fp;

    snprintf(path, cap, "/tmp/osh_scoring_save_detect_%ld.tmp", (long) getpid());
    fp = fopen(path, "w");
    ASSERT_TRUE(fp != NULL);
    ASSERT_TRUE(fputs(content, fp) >= 0);
    ASSERT_TRUE(fclose(fp) == 0);
}
