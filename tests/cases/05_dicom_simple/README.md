- TEST 05_DICOM_SIMPLE -

- DESCRIPTION:
Loads the DCPT_headphantom CT (512x512x177 voxels, 0.5859 mm pixel spacing,
1.5 mm slice spacing) from tests/fixtures/dicom/DCPT_headphantom.

The fixture is a real clinical proton plan (Brain_fin2) with three fields, all
at gantry=0 / couch=0.  The default geo.dat / beam.dat are a simplified
single-field approximation (126 MeV, square field).  Per-field geo.dat files
(geo_field1.dat .. geo_field3.dat) carry the correct RTPLAN-derived isocenter
placement; matching beam.dat files are TODO.

- RTPLAN: Brain_fin2 (RN.1.2.246.352.71.5...265919....dcm)

  Isocenter PCS (all fields): (0.0, -170.159, -2.122) mm
  RangeShifter (all fields):  WET 57.0 mm, ICRU 219 (polycarbonate/Lexan), Setting IN
                              IsocenterToRangeShifterDistance = distance to RS midplane.
                              Physical thickness 50 mm assumed (WET 57 mm at ~1.14 WER).

  Field 1
    GantryAngle: 0 deg   CouchAngle: 0 deg
    NominalEnergy: 186.197 MeV
    IsocenterToRangeShifter: 275.531 mm   SnoutPosition: 232.531 mm
    RS midplane z: +27.553 cm   RS slab z: [+25.053, +30.053] cm
    RTDOSE: RD.1.2.246.352.71.7.37402163639.2864312.20240227185700.dcm

  Field 2
    GantryAngle: 0 deg   CouchAngle: 0 deg
    NominalEnergy: 156.92 MeV
    IsocenterToRangeShifter: 293.617 mm   SnoutPosition: 250.617 mm
    RS midplane z: +29.362 cm   RS slab z: [+26.862, +31.862] cm
    RTDOSE: RD.1.2.246.352.71.7.37402163639.2864313.20240227185700.dcm

  Field 3
    GantryAngle: 0 deg   CouchAngle: 0 deg
    NominalEnergy: 154.114 MeV
    IsocenterToRangeShifter: 252.235 mm   SnoutPosition: 209.235 mm
    RS midplane z: +25.223 cm   RS slab z: [+22.723, +27.723] cm
    RTDOSE: RD.1.2.246.352.71.7.37402163639.2864314.20240227185700.dcm

- DCM PLACEMENT (geo.dat and geo_field*.dat):
  The simulation universe follows IEC 61217:
    Universe X = patient left
    Universe Y = cranial direction (head-first patients)
    Universe Z = patient anterior = nozzle direction at gantry 0

  The DCM card takes the isocenter in DICOM LPS patient coordinates [mm]:
    DCM CTBOX <ct_dir> HFS <gantry_deg> <couch_deg> <iso_x_mm> <iso_y_mm> <iso_z_mm>

  DICOM origin: (-149.707, -316.207, -139.0) mm  (LPS patient coords)
  Isocenter PCS (all fields): (0.0, -170.159, -2.122) mm

  geo.dat uses gantry=90 / couch=+30 as a non-trivial rotation smoke test:
    DCM CTBOX ../../fixtures/dicom/DCPT_headphantom HFS 90.0 30.0 \
        0.0 -170.15853658537 -2.1219512195122

  The app converts iso_mm to universe tx/ty/tz [cm] via the patient-position
  base rotation; callers no longer need to do this conversion manually.

- CT universe extent (HFS, all three fields identical):
    universe X: -14.971 to +15.029 cm   (512 x 0.5859 mm, along DICOM X = patient left-right)
    universe Y: -13.688 to +12.862 cm   (177 x 1.5 mm, along DICOM Z = cranial-caudal)
    universe Z: -15.395 to +14.605 cm   (512 x 0.5859 mm, along -DICOM Y = ant-post = beam dir)

- PURPOSE:
End-to-end validation of CT loading, HU-to-material mapping (Schneider2000),
voxel transport, Mesh scoring, and RTDOSE scoring/write-back.  Run with
--dry-run by default (validates parsing + compile only).

To enable transport, add args.cmake in this directory:
  set(OSH_ARGS "-v" "--outdir" "${WORK_DIR}")

- GEOMETRY:
  BHBOX : outer blackhole RPP  (universe X/Z: +-25 cm, Y: +-20 cm)
  OUTER : vacuum buffer RPP    (universe X/Z: +-20 cm, Y: -14.5..+13.5 cm)
  CTBOX : DCM head phantom     (HFS patient position, gantry=90, couch=+30)
  DCM card format: name ct_dir patient_pos gantry_deg couch_deg iso_x_mm iso_y_mm iso_z_mm
    iso_{x,y,z}_mm: treatment isocenter in DICOM LPS patient coordinates [mm]
    (app converts to universe tx/ty/tz internally via the patient-position base rotation)

- MATERIAL:
  HUTABLE Schneider2000 (24-bin CT calibration)
  CT zone assigned schneider_00 as placeholder (HU LUT overrides at runtime)

- BEAM:
  126 MeV protons, monoenergetic
  5x5 cm square field (BEAMSIGMA -2.5 -2.5)
  Direction: (0, 0, -1)  [BEAMDIR 180: nozzle at universe +Z, beam travels in -Z]
  BEAMPOS -15 cm (PZALIGN/beam-local frame) = universe z=+15 cm start position
  Universe +Z is the nozzle direction at gantry 0 per IEC 61217 (anterior face
  of the patient for HFS supine).  Beam enters the head phantom anteriorly.

- SCORING OUTPUTS:
  NB_ddc.dat  — depth-dose profile along beam axis (Mesh, 1x1x270 bins, TEXT)
  NB_XZ.bdo   — XZ plane fluence/energy (Mesh, 300x1x300, BDO)
  NB_ct.bdo   — energy scored onto the full CT voxel grid (DicomCT → Mesh,
                512x512x177 bins, BDO; ~355 MB)
                Current limitation: this path assumes zero CT rotation.
  NB_rtdose.dcm — energy scored onto the RTDOSE grid from fixture RD*.dcm
                  (DicomRTDOSE geometry → Mesh conversion at app level;
                   RTDOSE round-trip write: template metadata preserved,
                   pixel data replaced with scored values / nstat /
                   dose_grid_scaling)

  Note: NB_rtdose.dcm stores energy-per-primary in dose_grid_scaling units,
  not absorbed dose.  Proper DOSE scoring (with rho x ds weighting) is deferred.

- RTDOSE FIXTURE:
  RD.1.2.246.352.71.7.37402163639.2864312.20240227185700.dcm
  300 cols x 144 rows x 354 frames, 1x1 mm pixel, 0.75 mm frame spacing
  Origin: (-149.707, -261.793, -139.000) mm (DICOM patient coords)
  In simulation world coords (after patient->world offset): x -15..+15 cm,
  y -9.58..+4.82 cm, z -13.24..+13.31 cm — beam passes through this volume.

- PLOTTING:
  After transport, overlay scored outputs on the CT with plot_dicom.py.
  The tool auto-detects BDO and DICOM RTDOSE formats; pass the CT series
  directory (not a glob) as the first argument.

  XZ scoring plane (BDO):
    python3 tools/plot_dicom.py \
      tests/fixtures/dicom/DCPT_headphantom \
      tests/cases/05_dicom_simple/NB_XZ.bdo \
      --geo tests/cases/05_dicom_simple/geo.dat \
      --log \
      -o tests/cases/05_dicom_simple/NB_XZ_overlay.png

  RTDOSE output (central XZ slice through the RTDOSE grid):
    python3 tools/plot_dicom.py \
      tests/fixtures/dicom/DCPT_headphantom/ \
      tests/cases/05_dicom_simple/NB_rtdose.dcm

  Add -o <file.png> to save instead of displaying interactively.
  Add --log for log10 scale.  Use --plane and --cut to choose a
  different projection or slice position.

  Requires Python packages: numpy, matplotlib, pydicom.

- PROFILING (cache-miss comparison):
  beam_pz.dat is a +Z reference beam (no BEAMDIR, BEAMPOS 0 0 -14) for
  comparing transport speed and LLC miss rate against the default -Z beam.

  The -b/-g/-m/-d short options override individual input files:
    build/bin/openshieldhit -n 10000 tests/cases/05_dicom_simple/
    build/bin/openshieldhit -n 10000 -b tests/cases/05_dicom_simple/beam_pz.dat \
        tests/cases/05_dicom_simple/

  With perf (requires linux-perf package):
    perf stat -e LLC-loads,LLC-load-misses \
        build/bin/openshieldhit -n 10000 tests/cases/05_dicom_simple/
    perf stat -e LLC-loads,LLC-load-misses \
        build/bin/openshieldhit -n 10000 -b tests/cases/05_dicom_simple/beam_pz.dat \
        tests/cases/05_dicom_simple/

  Expected: -Z beam (default) shows higher LLC-load-misses than +Z beam due to
  backward traversal through the Z-major CT HU array (512x512x2 = 512 KB per
  Z-slice stride, defeating the hardware prefetcher).
