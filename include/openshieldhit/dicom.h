#ifndef OPENSHIELDHIT_DICOM_H
#define OPENSHIELDHIT_DICOM_H

#include <stddef.h>
#include <stdint.h>

#include "openshieldhit/diag.h"
#include "openshieldhit/status.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief CT image series read from a directory of per-slice DICOM files.
 *
 * @details
 * Slices are sorted by z-position. All coordinates are in mm, following the
 * DICOM Patient Coordinate System.
 */
struct osh_dicom_ct {
    double origin[3];        /**< DICOM Image Position Patient (center of first voxel in first slice) [mm] */
    double row_cosine[3];    /**< Image Orientation Patient: row direction */
    double col_cosine[3];    /**< Image Orientation Patient: column direction */
    double pixel_spacing[2]; /**< [row_spacing, col_spacing] in mm */
    double slice_spacing;    /**< z-spacing between slices [mm] */
    int rows;
    int cols;
    int n_slices;
    double rescale_intercept; /**< HU = pixel * rescale_slope + rescale_intercept */
    double rescale_slope;
    int16_t *pixels; /**< [n_slices * rows * cols], freed by osh_dicom_ct_free */
};

/**
 * @brief RTDOSE file: readable and round-trip pixel-writable.
 *
 * @details
 * @p pixels points directly into @p _raw. Modify pixels in-place, then call
 * osh_dicom_rtdose_write() to write the modified file back. Only pixel values
 * may be changed; array dimensions and all metadata are fixed.
 *
 * Dose in Gy: `dose = pixels[i] * dose_grid_scaling`.
 */
struct osh_dicom_rtdose {
    double origin[3];
    double row_cosine[3];
    double col_cosine[3];
    double pixel_spacing[2];
    double *frame_offsets; /**< z-offsets of dose planes [mm], n_frames entries */
    int rows;
    int cols;
    int n_frames;
    double dose_grid_scaling;
    uint32_t *pixels; /**< [n_frames * rows * cols], points into _raw */
    size_t n_pixels;
    /* private — round-trip support, do not modify directly */
    unsigned char *_raw;
    size_t _raw_size;
    size_t _pixel_data_offset;
};

/**
 * @brief Read a CT series from a directory of per-slice DICOM files.
 *
 * @param[in]  dir   Directory path containing the CT .dcm files.
 * @param[out] ct    Populated on success; free with osh_dicom_ct_free().
 * @param[in]  diag  Diagnostics sink; NULL means silent.
 */
enum osh_status osh_dicom_ct_read(char const *dir, struct osh_dicom_ct *ct, struct osh_diag_sink const *diag);

void osh_dicom_ct_free(struct osh_dicom_ct *ct);

/**
 * @brief Read a single RTDOSE file.
 *
 * @param[in]  path  Path to the RTDOSE .dcm file.
 * @param[out] rd    Populated on success; free with osh_dicom_rtdose_free().
 * @param[in]  diag  Diagnostics sink; NULL means silent.
 */
enum osh_status osh_dicom_rtdose_read(char const *path, struct osh_dicom_rtdose *rd, struct osh_diag_sink const *diag);

/**
 * @brief Write a modified RTDOSE file back to disk.
 *
 * @details
 * Only pixel values (already modified in-place via @p rd->pixels) are
 * written; the DICOM metadata is unchanged.
 *
 * @param[in] path  Output path (may be the same as the source file).
 * @param[in] rd    RTDOSE struct whose pixels have been modified.
 * @param[in] diag  Diagnostics sink; NULL means silent.
 */
enum osh_status
osh_dicom_rtdose_write(char const *path, struct osh_dicom_rtdose const *rd, struct osh_diag_sink const *diag);

void osh_dicom_rtdose_free(struct osh_dicom_rtdose *rd);

#ifdef __cplusplus
}
#endif

#endif /* OPENSHIELDHIT_DICOM_H */
