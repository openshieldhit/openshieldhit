#include "gemca/voxel/osh_gemca2_voxel_parse.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/osh_logger.h"
#include "common/osh_readline.h"
#include "gemca/voxel/osh_gemca2_voxel.h"
#include "gemca/voxel/osh_gemca2_voxel_defines.h"
#include "gemca/voxel/osh_gemca2_voxel_keys.h"

static int _parse_header(struct oshfile *shf, struct voxelct *ct);
static int _load_ctx(struct voxelct *ct);
static inline int16_t _swap_int16(int16_t val);

int osh_gemca_voxel_load(char const *fname, struct voxelct *ct) {

    char *fhed;
    char *fctx;
    int fhed_l;
    int fctx_l;

    struct shfile *shf;

    shf = calloc(1, sizeof(struct oshfile));
    if (shf == NULL) {
        osh_err(EX_UNAVAILABLE, "could not allocate memory.\n");
    }

    /* prepare the .hed filename */
    fhed_l = strlen(fname) + strlen(OSH_GEMCA_VOXEL_SUFFIX_HED);
    fctx_l = strlen(fname) + strlen(OSH_GEMCA_VOXEL_SUFFIX_CTX);

    fhed = calloc(fhed_l + 1, sizeof(char));
    fctx = calloc(fctx_l + 1, sizeof(char));

    snprintf(fhed, fhed_l, "%s%s", fname, OSH_GEMCA_VOXEL_SUFFIX_HED);
    snprintf(fctx, fctx_l, "%s%s", fname, OSH_GEMCA_VOXEL_SUFFIX_CTX);

    ct->fname_hed = fhed;
    ct->fname_ctx = fctx;

    shf = osh_fopen(fhed);
    _parse_header(shf, ct);
    osh_fclose(&shf);
    free(shf);

    _load_ctx(ct);

    return 1;
}

/* parse the .hed file */
/* ct->filename must have been set prior */
int _parse_header(struct oshfile *shf, struct voxelct *ct) {
    char *key = NULL;
    char *args = NULL;
    char *line = NULL;

    int j;
    unsigned int ui;
    double a, b, c;

    int lineno;
    int len;

    while (osh_readline_key(shf, &line, &key, &args, &lineno) > 0) { /* version */
        len = strlen(args);
        if (strcasecmp(OSH_GEMCA_VOXEL_KEY_VERSION, key) == 0) {
            ct->version = calloc(len + 1, sizeof(char));
            strncpy(ct->version, args, len + 1);
        }

        else if (strcasecmp(OSH_GEMCA_VOXEL_KEY_MODALITY, key) == 0) { /* modality */
            ct->modality = calloc(len + 1, sizeof(char));
            strncpy(ct->modality, args, len + 1);
        }

        else if (strcasecmp(OSH_GEMCA_VOXEL_KEY_CREATEDBY, key) == 0) { /* created_by */
            ct->created_by = calloc(len + 1, sizeof(char));
            strncpy(ct->created_by, args, len + 1);
        }

        else if (strcasecmp(OSH_GEMCA_VOXEL_KEY_CINFO, key) == 0) { /* creation_info */
            ct->creation_info = calloc(len + 1, sizeof(char));
            strncpy(ct->creation_info, args, len + 1);
        }

        else if (strcasecmp(OSH_GEMCA_VOXEL_KEY_PRIMARYVIEW, key) == 0) { /* primary_view */
            ct->primary_view = calloc(len + 1, sizeof(char));
            strncpy(ct->primary_view, args, len + 1);
        }

        else if (strcasecmp(OSH_GEMCA_VOXEL_KEY_DATATYPE, key) == 0) { /* data_type */
            if (strcasecmp(OSH_GEMCA_VOXEL_KEY_INTEGER, args) == 0) {
                ct->data_type = OSH_GEMCA_VOXEL_INTEGER;
            } else if (strcasecmp(OSH_GEMCA_VOXEL_KEY_FLOAT, args) == 0) {
                ct->data_type = OSH_GEMCA_VOXEL_FLOAT;
            } else {
                osh_err(EX_CONFIG,
                        "_parse_header(): unknown data_type '%s' in file %s line %lu.",
                        args,
                        ct->fname_hed,
                        lineno);
            }
        }

        else if (strcasecmp(OSH_GEMCA_VOXEL_KEY_NUMBYTES, key) == 0) { /* num_bytes */
            ct->data_type_size = atoi(args);
        }

        else if (strcasecmp(OSH_GEMCA_VOXEL_KEY_BYTEORDER, key) == 0) { /* byte_order */
            if (strcasecmp(OSH_GEMCA_VOXEL_KEY_AIX, args) == 0) {
                ct->byte_order = OSH_GEMCA_VOXEL_BIGENDIAN;
            } else if (strcasecmp(OSH_GEMCA_VOXEL_KEY_VMS, args) == 0) {
                ct->byte_order = OSH_GEMCA_VOXEL_LITTLEENDIAN;
            } else {
                osh_err(EX_CONFIG,
                        "_parse_header(): unknown byte_order '%s' in file %s line %lu.",
                        args,
                        ct->fname_hed,
                        lineno);
            }
        }

        else if (strcasecmp(OSH_GEMCA_VOXEL_KEY_PATIENTNAME, key) == 0) { /* patient_name */
            ct->patient_name = calloc(len + 1, sizeof(char));
            strncpy(ct->patient_name, args, len + 1);
        }

        else if (strcasecmp(OSH_GEMCA_VOXEL_KEY_SLICEDIM, key) == 0) { /* slice_dimension */
            ct->slice_dimension = atoi(args);
        }

        else if (strcasecmp(OSH_GEMCA_VOXEL_KEY_PIXELSIZE, key) == 0) {   /* pixel_size */
            ct->pixel_size = atof(args) * 0.1;                            /* [mm] -> [cm] */
        } else if (strcasecmp(OSH_GEMCA_VOXEL_KEY_SLICEDIST, key) == 0) { /* slice_distance */
            ct->slice_distance = atof(args) * 0.1;                        /* [mm] -> [cm] */
        }

        else if (strcasecmp(OSH_GEMCA_VOXEL_KEY_SLICENUM, key) == 0) { /* slice_number (= number fo slices...) */
            ct->slice_number = atoi(args);
        }

        else if (strcasecmp(OSH_GEMCA_VOXEL_KEY_XOFFSET, key) == 0) { /* OFFSET */
            ct->offset[0] = atoi(args);
        } else if (strcasecmp(OSH_GEMCA_VOXEL_KEY_YOFFSET, key) == 0) {
            ct->offset[1] = atoi(args);
        } else if (strcasecmp(OSH_GEMCA_VOXEL_KEY_ZOFFSET, key) == 0) {
            ct->offset[2] = atoi(args);
        }

        else if (strcasecmp(OSH_GEMCA_VOXEL_KEY_DIMX, key) == 0) { /* DIMx,y,z */
            ct->dim[0] = atoi(args);
        } else if (strcasecmp(OSH_GEMCA_VOXEL_KEY_DIMY, key) == 0) {
            ct->dim[1] = atoi(args);
        } else if (strcasecmp(OSH_GEMCA_VOXEL_KEY_DIMZ, key) == 0) {
            ct->dim[2] = atoi(args);
        }

        else if (strcasecmp(OSH_GEMCA_VOXEL_KEY_ZTABLE, key) == 0) { /* z_table */
            ct->has_ztable = 0;
            if (strcasecmp(OSH_GEMCA_VOXEL_KEY_YES, args) == 0) {
                ct->has_ztable = 1;
                /* allocate memory for the z_table */
                ct->ztable_pos = calloc(ct->slice_number, sizeof(double)); // TODO: ask why not DIMZ
                ct->ztable_thickness = calloc(ct->slice_number, sizeof(double));
                ct->ztable_gantry_tilt = calloc(ct->slice_number, sizeof(double));
            } else if (strcasecmp(OSH_GEMCA_VOXEL_KEY_NO, args) == 0) {
                ct->has_ztable = 0;
            } else {
                osh_err(EX_CONFIG,
                        "_parse_header(): did not understand z_table '%s' in file %s line %lu.",
                        args,
                        ct->fname_hed,
                        lineno);
            }
        } else if (strcasecmp(OSH_GEMCA_VOXEL_KEY_SLICENO, key) == 0) { /* slice_no: z_table data */
            /* we now need to read a z-table */
            if (!ct->has_ztable) {
                osh_warn("_parse_header(): in %s line %lu. No 'z_table yes', skipping the rest of the file. ",
                         ct->fname_hed,
                         lineno);
                break;
            }
            while (osh_readline(shf, &line, &lineno) > 0) { /* read the z_table */
                j = sscanf(line, "%ui %lf %lf %lf\n", &ui, &a, &b, &c);

                if (j != 4) { /* check if the number of columns is right in the z_table */
                    osh_err(EX_CONFIG,
                            "_parse_header(): z_table wrong number of columns (should be 4) in file %s line %lu.",
                            ct->fname_hed,
                            lineno);
                }

                if (ui > ct->slice_number) { /* check bounds */
                    osh_err(EX_CONFIG,
                            "_parse_header(): z_table slice index larger than slice_number in file %s line %lu.",
                            ct->fname_hed,
                            lineno);
                }

                if (ui < 1) { /* check bounds */
                    osh_err(EX_CONFIG,
                            "_parse_header(): z_table slice index must start at 1. In file %s line %lu.",
                            ct->fname_hed,
                            lineno);
                }
                --ui; /* z_table index starts at 1. Our arrays start at 0 */
                ct->ztable_pos[ui] = a;
                ct->ztable_thickness[ui] = b;
                ct->ztable_gantry_tilt[ui] = c;
            } /* end of z-table read while loop */
            free(line);
        } /* end of SLICENO else if */

        free(line);
    } /* end of while loop */
    free(line);
    return 1;
}

/* parse the .hed file */
static int _load_ctx(struct voxelct *ct) {
    size_t n; /* number of elements in CT cube */
    FILE *fp;

    size_t i;

    n = ct->dim[0] * ct->dim[1] * ct->dim[2];
    // ct->data_type_size

    /* allocate memory for CT data */

    if (ct->data_type_size != 2) {
        osh_err(EX_CONFIG, "CTX data_type must be integer 2 byte (signed short)");
    }

    if (ct->data_type == OSH_GEMCA_VOXEL_INTEGER)
        ct->hu = calloc(ct->dim[0] * ct->dim[1] * ct->dim[2], sizeof(ct->hu));

    fp = fopen(ct->fname_hed, "rb");
    while (!feof(fp)) {
        fread(ct->hu, sizeof(ct->hu), n, fp);
    }
    fclose(fp);

    /* check if we need to bytewap */
    if (ct->byte_order == OSH_GEMCA_VOXEL_BIGENDIAN) {
        for (i = 0; i < n; i++) {
            _swap_int16(ct->hu[i]);
        }
    }

    return 1;
}

// Byte swap short
// https://stackoverflow.com/questions/2182002/convert-big-endian-to-little-endian-in-c-without-using-provided-func
static inline int16_t _swap_int16(int16_t val) {
    return (val << 8) | ((val >> 8) & 0xFF);
}
