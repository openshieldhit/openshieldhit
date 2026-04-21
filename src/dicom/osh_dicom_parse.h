#ifndef OSH_DICOM_PARSE_H
#define OSH_DICOM_PARSE_H

#include <stddef.h>
#include <stdint.h>

#include "openshieldhit/diag.h"

/**
 * @brief Tag callback for osh_dicom_walk().
 *
 * @param[in] group    Tag group number.
 * @param[in] element  Tag element number.
 * @param[in] vr       Two-character VR (not null-terminated).
 * @param[in] value    Pointer to the tag's value bytes within the source buffer.
 * @param[in] length   Value length in bytes.
 * @param[in] user     Caller context.
 * @return 1 to continue walking, 0 to stop early.
 */
typedef int (*osh_dicom_tag_fn)(
    uint16_t group, uint16_t element, char const vr[2], unsigned char const *value, uint32_t length, void *user);

/**
 * @brief Load a file into a heap buffer.
 * @return Allocated buffer (caller frees), or NULL on failure.
 */
unsigned char *osh_dicom_load_file(char const *path, size_t *out_size, struct osh_diag_sink const *diag);

/**
 * @brief Walk all top-level tags in a DICOM Part-10 buffer (explicit VR LE).
 *
 * @details
 * Does not recurse into sequences. Sufficient for all flat metadata tags
 * and pixel data at the top level of CT, RTDOSE, and RTSTRUCT files.
 *
 * @return 1 on success, 0 if the buffer is not a valid DICOM file.
 */
int osh_dicom_walk(
    unsigned char const *buf, size_t size, osh_dicom_tag_fn fn, void *user, struct osh_diag_sink const *diag);

/* Value decoders ---------------------------------------------------------- */

/** Parse a backslash-separated DS (Decimal String) array. Returns count parsed. */
int osh_dicom_ds_array(unsigned char const *value, uint32_t length, double *out, int max_count);

/** Parse a single DS value. */
double osh_dicom_ds(unsigned char const *value, uint32_t length);

/** Parse a US (Unsigned Short) value. */
uint16_t osh_dicom_us(unsigned char const *value, uint32_t length);

#endif /* OSH_DICOM_PARSE_H */
