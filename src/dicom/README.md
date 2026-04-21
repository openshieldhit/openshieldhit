# DICOM

This directory contains a deliberately small DICOM reader/writer layer for the
subset of medical-image data that OpenShieldHIT currently cares about.

## Scope

- CT image series: read a directory of per-slice DICOM files into one dense
  volume.
- RTDOSE: read one dose file, expose its geometry and pixel grid, and write the
  modified pixel payload back to disk.

The goal is not to provide a general-purpose DICOM toolkit. The code here is a
small internal utility layer that supports inspection, import, and simple
round-trip dose editing.

## Current assumptions

- DICOM Part-10 file layout with `DICM` preamble
- Explicit-VR little-endian top-level tag walking
- No recursive sequence parsing
- Flat metadata extraction for CT and RTDOSE only

Those assumptions are enough for the current CT/RTDOSE workflows, but they are
not sufficient for arbitrary DICOM objects.

## Layout

- `osh_dicom_parse.c`
  Minimal file loader and top-level tag walker plus a few value decoders.

- `osh_dicom_ct.c`
  CT-series loader. Reads per-slice files, filters to CT modality, sorts by
  z-position, and assembles one dense voxel array.

- `osh_dicom_rtdose.c`
  RTDOSE reader/writer. Keeps the original raw file buffer so the pixel payload
  can be modified in place and written back without rebuilding the object.

## Public API

The installed public entry points are declared in
`include/openshieldhit/dicom.h`:

- `osh_dicom_ct_read()` / `osh_dicom_ct_free()`
- `osh_dicom_rtdose_read()` / `osh_dicom_rtdose_write()` /
  `osh_dicom_rtdose_free()`

All diagnostics are emitted through a caller-owned `struct osh_diag_sink`.
