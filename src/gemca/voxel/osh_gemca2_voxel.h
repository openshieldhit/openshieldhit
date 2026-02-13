#ifndef _OSH_GEMCA2_VOXEL
#define _OSH_GEMCA2_VOXEL

#include <stdint.h>

struct voxelct {
    char *fname_hed; /* filename */
    char *fname_ctx; /* filename */
    char *version;
    char *modality;
    char *created_by;    /* */
    char *creation_info; /* */
    char *primary_view;  /* */
    char *patient_name;  /* */
    char data_type;
    char data_type_size; /* data type size in bytes, aka the "num_bytes" string */
    char byte_order;     /* little or big endian */
    unsigned int slice_dimension;
    double pixel_size;     /* in cm */
    double slice_distance; /* in cm */
    unsigned int slice_number;
    int offset[3]; /* x-, y-, zoffset */
    unsigned int dim[3];
    char has_ztable;
    double *ztable_pos;
    double *ztable_thickness;
    double *ztable_gantry_tilt;
    double *data; /* CT data stored as float */
    int16_t *hu;
};

#endif /* _OSH_GEMCA2_VOXEL */
