#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "apps/osh/osh_scoring_parse_internal.h"
#include "common/osh_diag.h"
#include "openshieldhit/scoring.h"
#include "openshieldhit/status.h"
#include "test_assert.h"

static void test_name(void);
static void test_axis(void);
static void test_rotation(void);
static void test_rotation_identity(void);
static void test_zones(void);
static void test_inputpath(void);
static void test_body(void);
static void test_unknown_key_is_ok(void);
static void test_error_missing_args(void);

/* Split a space-separated string into words in-place.  Returns word count. */
static int tokenize(char *buf, char **words, int max_words) {
    int n = 0;
    char *p = buf;
    while (*p && n < max_words) {
        while (*p == ' ')
            p++;
        if (!*p)
            break;
        words[n++] = p;
        while (*p && *p != ' ')
            p++;
        if (*p)
            *p++ = '\0';
    }
    return n;
}

static struct osh_scoring_geometry_def make_geo(void) {
    struct osh_scoring_geometry_def geo;
    memset(&geo, 0, sizeof(geo));
    return geo;
}

static void free_geo(struct osh_scoring_geometry_def *geo) {
    size_t z; /* Zone selector index for freeing zone_names[]. */
    free(geo->kind);
    free(geo->name);
    free(geo->axes);
    for (z = 0u; z < geo->nzone_indices; ++z) {
        free(geo->zone_names[z]);
    }
    free((void *) geo->zone_names);
    free(geo->zone_indices);
    free(geo->zone_volumes);
    free(geo->vox_rtdose_path);
    free(geo->vox_body_name);
    geo->naxes = 0;
    geo->nzone_indices = 0u;
}

static int nearly(double a, double b) {
    return fabs(a - b) < 1e-9;
}

int main(void) {
    test_name();
    test_axis();
    test_rotation();
    test_rotation_identity();
    test_zones();
    test_inputpath();
    test_body();
    test_unknown_key_is_ok();
    test_error_missing_args();
    return 0;
}

static void test_name(void) {
    struct osh_scoring_geometry_def geo = make_geo();
    char line[] = "name my_detector";
    char *words[8];
    int nw = tokenize(line, words, 8);
    int found = 0;
    enum osh_status rc = osh_scoring_parse_geometry_line(&geo, NULL, words, nw, "test", 1, &found);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(found == 1);
    ASSERT_TRUE(geo.name != NULL);
    ASSERT_TRUE(strcmp(geo.name, "my_detector") == 0);
    free_geo(&geo);
}

static void test_axis(void) {
    struct osh_scoring_geometry_def geo = make_geo();
    char lx[] = "x -10.0 10.0 100";
    char ly[] = "y -5.0 5.0 50";
    char lz[] = "z 0.0 20.0 200";
    char *words[8];
    int nw;
    int found;

    nw = tokenize(lx, words, 8);
    ASSERT_TRUE(osh_scoring_parse_geometry_line(&geo, NULL, words, nw, "t", 1, &found) == OSH_OK);
    ASSERT_TRUE(found == 1);

    nw = tokenize(ly, words, 8);
    ASSERT_TRUE(osh_scoring_parse_geometry_line(&geo, NULL, words, nw, "t", 2, &found) == OSH_OK);

    nw = tokenize(lz, words, 8);
    ASSERT_TRUE(osh_scoring_parse_geometry_line(&geo, NULL, words, nw, "t", 3, &found) == OSH_OK);

    ASSERT_TRUE(geo.naxes == 3u);
    ASSERT_TRUE(strcmp(geo.axes[0].label, "X") == 0);
    ASSERT_TRUE(nearly(geo.axes[0].lo, -10.0) && nearly(geo.axes[0].hi, 10.0));
    ASSERT_TRUE(geo.axes[0].nbins == 100);
    ASSERT_TRUE(strcmp(geo.axes[1].label, "Y") == 0);
    ASSERT_TRUE(nearly(geo.axes[1].lo, -5.0) && nearly(geo.axes[1].hi, 5.0));
    ASSERT_TRUE(geo.axes[1].nbins == 50);
    ASSERT_TRUE(strcmp(geo.axes[2].label, "Z") == 0);
    ASSERT_TRUE(nearly(geo.axes[2].lo, 0.0) && nearly(geo.axes[2].hi, 20.0));
    ASSERT_TRUE(geo.axes[2].nbins == 200);
    free_geo(&geo);
}

static void test_rotation(void) {
    struct osh_scoring_geometry_def geo = make_geo();
    char line[] = "rotation 90.0 45.0";
    char *words[8];
    int nw = tokenize(line, words, 8);
    int found = 0;
    enum osh_status rc = osh_scoring_parse_geometry_line(&geo, NULL, words, nw, "test", 1, &found);
    ASSERT_TRUE(rc == OSH_OK);
    ASSERT_TRUE(found == 1);
    ASSERT_TRUE(geo.has_rotation == 1u);
    /* Verify rows of the 3x3 rotation block are unit vectors */
    {
        int i;
        for (i = 0; i < 3; i++) {
            double n2 = geo.t[i * 4 + 0] * geo.t[i * 4 + 0] + geo.t[i * 4 + 1] * geo.t[i * 4 + 1]
                        + geo.t[i * 4 + 2] * geo.t[i * 4 + 2];
            ASSERT_TRUE(nearly(n2, 1.0));
        }
    }
    /* Translation column must be zero */
    ASSERT_TRUE(nearly(geo.t[3], 0.0));
    ASSERT_TRUE(nearly(geo.t[7], 0.0));
    ASSERT_TRUE(nearly(geo.t[11], 0.0));
    free_geo(&geo);
}

static void test_rotation_identity(void) {
    struct osh_scoring_geometry_def geo = make_geo();
    char line[] = "rot 0.0 0.0";
    char *words[8];
    int nw = tokenize(line, words, 8);
    ASSERT_TRUE(osh_scoring_parse_geometry_line(&geo, NULL, words, nw, "test", 1, NULL) == OSH_OK);
    ASSERT_TRUE(geo.has_rotation == 1u);
    /* Identity: diagonal ones, off-diagonal zeros */
    ASSERT_TRUE(nearly(geo.t[0], 1.0) && nearly(geo.t[5], 1.0) && nearly(geo.t[10], 1.0));
    ASSERT_TRUE(nearly(geo.t[1], 0.0) && nearly(geo.t[2], 0.0));
    ASSERT_TRUE(nearly(geo.t[4], 0.0) && nearly(geo.t[6], 0.0));
    ASSERT_TRUE(nearly(geo.t[8], 0.0) && nearly(geo.t[9], 0.0));
    free_geo(&geo);
}

static void test_zones(void) {
    struct osh_scoring_geometry_def geo = make_geo();
    char line0[] = "zone WaterBox";
    char line0v[] = "volume 13.37";
    char line1[] = "zone Target";
    char *words[8];
    int nw = tokenize(line0, words, 8);
    int found = 0;
    ASSERT_TRUE(osh_scoring_parse_geometry_line(&geo, NULL, words, nw, "test", 1, &found) == OSH_OK);
    ASSERT_TRUE(found == 1);

    nw = tokenize(line0v, words, 8);
    ASSERT_TRUE(osh_scoring_parse_geometry_line(&geo, NULL, words, nw, "test", 2, &found) == OSH_OK);

    nw = tokenize(line1, words, 8);
    ASSERT_TRUE(osh_scoring_parse_geometry_line(&geo, NULL, words, nw, "test", 3, &found) == OSH_OK);

    /* The parser stores zone selectors by name verbatim; resolution to transport
     * indices happens later at app level (osh_scoring_resolve_zone_names). */
    ASSERT_TRUE(geo.nzone_indices == 2u);
    ASSERT_TRUE(strcmp(geo.zone_names[0], "WaterBox") == 0);
    ASSERT_TRUE(strcmp(geo.zone_names[1], "Target") == 0);
    ASSERT_TRUE(nearly(geo.zone_volumes[0], 13.37));
    ASSERT_TRUE(nearly(geo.zone_volumes[1], 0.0));
    free_geo(&geo);
}

static void test_inputpath(void) {
    struct osh_scoring_geometry_def geo = make_geo();
    char line[] = "inputpath /some/path/dose.dcm";
    char *words[8];
    int nw = tokenize(line, words, 8);
    int found = 0;
    ASSERT_TRUE(osh_scoring_parse_geometry_line(&geo, NULL, words, nw, "test", 1, &found) == OSH_OK);
    ASSERT_TRUE(found == 1);
    ASSERT_TRUE(geo.vox_rtdose_path != NULL);
    ASSERT_TRUE(strcmp(geo.vox_rtdose_path, "/some/path/dose.dcm") == 0);
    free_geo(&geo);
}

static void test_body(void) {
    struct osh_scoring_geometry_def geo = make_geo();
    char line[] = "body CTBOX";
    char *words[8];
    int nw = tokenize(line, words, 8);
    int found = 0;
    ASSERT_TRUE(osh_scoring_parse_geometry_line(&geo, NULL, words, nw, "test", 1, &found) == OSH_OK);
    ASSERT_TRUE(found == 1);
    ASSERT_TRUE(geo.vox_body_name != NULL);
    ASSERT_TRUE(strcmp(geo.vox_body_name, "CTBOX") == 0);
    free_geo(&geo);
}

static void test_unknown_key_is_ok(void) {
    struct osh_scoring_geometry_def geo = make_geo();
    char line[] = "SomeUnknownKey arg1";
    char *words[8];
    int nw = tokenize(line, words, 8);
    int found = 0;
    /* Unknown keys must return OSH_OK with found=0 (dispatcher ignores them). */
    ASSERT_TRUE(osh_scoring_parse_geometry_line(&geo, NULL, words, nw, "test", 1, &found) == OSH_OK);
    ASSERT_TRUE(found == 0);
    free_geo(&geo);
}

static void test_error_missing_args(void) {
    struct osh_scoring_geometry_def geo;
    char *words[8];
    int nw;
    enum osh_status rc;

    /* name with no argument */
    memset(&geo, 0, sizeof(geo));
    {
        char l[] = "name";
        nw = tokenize(l, words, 8);
    }
    rc = osh_scoring_parse_geometry_line(&geo, NULL, words, nw, "test", 1, NULL);
    ASSERT_TRUE(rc == OSH_EPARSE);

    /* axis with too few arguments */
    memset(&geo, 0, sizeof(geo));
    {
        char l[] = "x -10.0 10.0";
        nw = tokenize(l, words, 8);
    }
    rc = osh_scoring_parse_geometry_line(&geo, NULL, words, nw, "test", 1, NULL);
    ASSERT_TRUE(rc == OSH_EPARSE);

    /* rotation with no arguments */
    memset(&geo, 0, sizeof(geo));
    {
        char l[] = "rotation";
        nw = tokenize(l, words, 8);
    }
    rc = osh_scoring_parse_geometry_line(&geo, NULL, words, nw, "test", 1, NULL);
    ASSERT_TRUE(rc == OSH_EPARSE);

    /* zone with no argument */
    memset(&geo, 0, sizeof(geo));
    {
        char l[] = "zone";
        nw = tokenize(l, words, 8);
    }
    rc = osh_scoring_parse_geometry_line(&geo, NULL, words, nw, "test", 1, NULL);
    ASSERT_TRUE(rc == OSH_EPARSE);

    /* inputpath with no argument */
    memset(&geo, 0, sizeof(geo));
    {
        char l[] = "inputpath";
        nw = tokenize(l, words, 8);
    }
    rc = osh_scoring_parse_geometry_line(&geo, NULL, words, nw, "test", 1, NULL);
    ASSERT_TRUE(rc == OSH_EPARSE);

    /* body with no argument */
    memset(&geo, 0, sizeof(geo));
    {
        char l[] = "body";
        nw = tokenize(l, words, 8);
    }
    rc = osh_scoring_parse_geometry_line(&geo, NULL, words, nw, "test", 1, NULL);
    ASSERT_TRUE(rc == OSH_EPARSE);
}
