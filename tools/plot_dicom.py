#!/usr/bin/env python3
"""Overlay an OpenShieldHIT scoring output on a DICOM CT slice.

Accepted scoring formats
------------------------
* OpenShieldHIT BDO2019 files  (*.bdo)
* DICOM RTDOSE files            (*.dcm)

The CT input may be a single DICOM file or a directory; when a directory is
given all files with PixelData are loaded as a series.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import struct
from dataclasses import dataclass
from pathlib import Path
from typing import Any, TypeAlias, cast

os.environ.setdefault("MPLCONFIGDIR", "/tmp/openshieldhit-matplotlib")

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt


OSHBDO_GEO_TYPE = 0xE000
OSHBDO_GEO_NAME = 0xE001
OSHBDO_GEO_P = 0xE002
OSHBDO_GEO_Q = 0xE003
OSHBDO_GEO_N = 0xE004
OSHBDO_PAG_TYPE = 0xDD30
OSHBDO_PAG_DATA = 0xDDBB
OSHBDO_PAG_DATA_UNIT = 0xDDBC

FloatArray: TypeAlias = npt.NDArray[np.float64]
IntArray: TypeAlias = npt.NDArray[np.int_]


@dataclass
class BdoData:
    path: Path
    geo_type: str
    geo_name: str
    p: FloatArray        # [x_lo, y_lo, z_lo] cm
    q: FloatArray        # [x_hi, y_hi, z_hi] cm
    n: IntArray          # [nx, ny, nz]
    pages: list[FloatArray]
    page_types: list[int]
    page_units: list[str]


@dataclass
class CtVolume:
    data_hu: FloatArray
    origin_cm: FloatArray
    spacing_cm: FloatArray


def _decode_c_string(payload: bytes) -> str:
    return payload.split(b"\0", 1)[0].decode("utf-8", errors="replace")


def _payload_dtype(pltype: str) -> np.dtype[np.generic] | None:
    pltype = pltype.strip("\0")
    if not pltype or pltype.startswith("S"):
        return None
    return np.dtype(pltype)


def read_bdo2019(path: Path) -> BdoData:
    tokens: dict[int, list[Any]] = {}
    endian = "<"

    with path.open("rb") as fp:
        preamble = fp.read(24)
        if len(preamble) != 24 or preamble[:6] != b"xSH12A":
            raise ValueError(f"{path}: not an OpenShieldHIT/SH12A BDO2019 file")
        if preamble[6:8] == b"MM":
            endian = ">"
        elif preamble[6:8] != b"II":
            raise ValueError(f"{path}: unsupported BDO endian marker {preamble[6:8]!r}")

        while True:
            header = fp.read(24)
            if not header:
                break
            if len(header) != 24:
                raise ValueError(f"{path}: truncated BDO token header")

            tag, raw_pltype, length = struct.unpack(f"{endian}Q8sQ", header)
            pltype = raw_pltype.split(b"\0", 1)[0].decode("ascii", errors="replace")
            if pltype.startswith("S"):
                payload_len = int(pltype[1:]) * int(length)
                value = _decode_c_string(fp.read(payload_len))
            else:
                dtype = _payload_dtype(pltype)
                if dtype is None:
                    payload_len = 0
                    value = None
                else:
                    payload_len = dtype.itemsize * int(length)
                    payload = fp.read(payload_len)
                    if len(payload) != payload_len:
                        raise ValueError(f"{path}: truncated BDO token payload")
                    value = np.frombuffer(payload, dtype=dtype).copy()
            tokens.setdefault(tag, []).append(value)

    try:
        p = cast(FloatArray, np.asarray(tokens[OSHBDO_GEO_P][-1], dtype=np.float64))
        q = cast(FloatArray, np.asarray(tokens[OSHBDO_GEO_Q][-1], dtype=np.float64))
        n = cast(IntArray, np.asarray(tokens[OSHBDO_GEO_N][-1], dtype=np.int_))
    except KeyError as exc:
        raise ValueError(f"{path}: missing required BDO geometry token") from exc

    pages = [cast(FloatArray, np.asarray(page, dtype=np.float64)) for page in tokens.get(OSHBDO_PAG_DATA, [])]
    if not pages:
        raise ValueError(f"{path}: no page data found")

    return BdoData(
        path=path,
        geo_type=str(tokens.get(OSHBDO_GEO_TYPE, [""])[-1]),
        geo_name=str(tokens.get(OSHBDO_GEO_NAME, [""])[-1]),
        p=p,
        q=q,
        n=n,
        pages=pages,
        page_types=[int(v[0]) for v in tokens.get(OSHBDO_PAG_TYPE, []) if v is not None and len(v)],
        page_units=[str(v) for v in tokens.get(OSHBDO_PAG_DATA_UNIT, [])],
    )


def read_rtdose_dicom(path: Path) -> BdoData:
    pydicom = _import_pydicom()
    ds = pydicom.dcmread(str(path))

    rows = int(ds.Rows)
    cols = int(ds.Columns)
    n_frames = int(ds.NumberOfFrames)

    col_spacing_cm = float(ds.PixelSpacing[1]) / 10.0
    row_spacing_cm = float(ds.PixelSpacing[0]) / 10.0
    frame_offsets = [float(v) for v in ds.GridFrameOffsetVector]
    z_spacing_cm = ((frame_offsets[1] - frame_offsets[0]) / 10.0) if len(frame_offsets) > 1 else col_spacing_cm

    ipp_cm = [float(v) / 10.0 for v in ds.ImagePositionPatient]
    lo_x = ipp_cm[0] - 0.5 * col_spacing_cm
    lo_y = ipp_cm[1] - 0.5 * row_spacing_cm
    lo_z = ipp_cm[2] + frame_offsets[0] / 10.0 - 0.5 * z_spacing_cm

    scaling = float(getattr(ds, "DoseGridScaling", 1.0))
    pixels = ds.pixel_array.astype(np.float64) * scaling  # (n_frames, rows, cols)

    return BdoData(
        path=path,
        geo_type="mesh",
        geo_name=str(getattr(ds, "SeriesDescription", path.stem)),
        p=cast(FloatArray, np.array([lo_x, lo_y, lo_z])),
        q=cast(FloatArray, np.array([lo_x + cols * col_spacing_cm,
                                     lo_y + rows * row_spacing_cm,
                                     lo_z + n_frames * z_spacing_cm])),
        n=cast(IntArray, np.array([cols, rows, n_frames])),
        pages=[pixels.ravel()],
        page_types=[0],
        page_units=["Gy"],
    )


def read_scoring(path: Path) -> BdoData:
    """Auto-detect BDO or DICOM RTDOSE and return a unified BdoData."""
    if path.suffix.lower() == ".dcm":
        return read_rtdose_dicom(path)
    try:
        return read_bdo2019(path)
    except ValueError:
        return read_rtdose_dicom(path)


def _import_pydicom() -> Any:
    try:
        import pydicom
    except ModuleNotFoundError as exc:
        raise SystemExit("error: pydicom is required; install it with `python3 -m pip install pydicom`") from exc
    return pydicom


def _slice_position(ds: Any) -> float:
    image_position = getattr(ds, "ImagePositionPatient", None)
    if image_position is not None:
        return float(image_position[2])
    return float(getattr(ds, "SliceLocation", getattr(ds, "InstanceNumber", 0)))


def read_ct_series(path: Path, origin_override_cm: FloatArray | None = None) -> CtVolume:
    pydicom = _import_pydicom()
    paths = [path] if path.is_file() else sorted(p for p in path.iterdir() if p.is_file())
    slices: list[Any] = []
    for candidate in paths:
        try:
            ds = pydicom.dcmread(str(candidate))
        except Exception:
            continue
        if hasattr(ds, "PixelData"):
            slices.append(ds)
    if not slices:
        raise ValueError(f"{path}: no DICOM image slices with PixelData found")

    ct_slices = [ds for ds in slices if str(getattr(ds, "Modality", "")).upper() == "CT"]
    if ct_slices:
        slices = ct_slices
    else:
        groups: dict[tuple[int, int, str], list[Any]] = {}
        for ds in slices:
            key = (
                int(getattr(ds, "Rows", -1)),
                int(getattr(ds, "Columns", -1)),
                str(getattr(ds, "SeriesInstanceUID", "")),
            )
            groups.setdefault(key, []).append(ds)
        slices = max(groups.values(), key=len)

    slices.sort(key=_slice_position)
    first = slices[0]
    spacing_mm = [float(first.PixelSpacing[1]), float(first.PixelSpacing[0])]
    if len(slices) > 1:
        z_positions = np.asarray([_slice_position(ds) for ds in slices], dtype=float)
        slice_spacing_mm = float(np.median(np.diff(z_positions)))
    else:
        slice_spacing_mm = float(getattr(first, "SliceThickness", 1.0))

    volume: list[FloatArray] = []
    for ds in slices:
        arr = np.asarray(ds.pixel_array, dtype=np.float64)
        slope = float(getattr(ds, "RescaleSlope", 1.0))
        intercept = float(getattr(ds, "RescaleIntercept", 0.0))
        volume.append(arr * slope + intercept)

    data_hu = cast(FloatArray, np.stack(volume, axis=0))
    spacing_cm = cast(
        FloatArray,
        np.asarray([0.1 * spacing_mm[0], 0.1 * spacing_mm[1], 0.1 * abs(slice_spacing_mm)], dtype=np.float64),
    )

    if origin_override_cm is not None:
        origin_cm = cast(FloatArray, origin_override_cm.astype(np.float64))
    elif getattr(first, "ImagePositionPatient", None) is not None:
        ipp_cm = cast(FloatArray, 0.1 * np.asarray(first.ImagePositionPatient, dtype=np.float64))
        origin_cm = ipp_cm - 0.5 * spacing_cm
    else:
        shape_cm = cast(FloatArray, np.asarray([data_hu.shape[2], data_hu.shape[1], data_hu.shape[0]], dtype=np.float64))
        origin_cm = -0.5 * spacing_cm * shape_cm

    return CtVolume(data_hu=data_hu, origin_cm=origin_cm, spacing_cm=spacing_cm)


# IEC 61217 patient-position base rotation matrices (rows = DICOM LPS axes in universe coords).
# Matches osh_patient_position_base_rotation() in src/common/osh_patient_position.c.
_PP_BASE_ROTATION: dict[str, FloatArray] = {
    "HFS":  cast(FloatArray, np.array([[ 1, 0, 0], [ 0, 0,-1], [ 0, 1, 0]], dtype=float)),
    "HFP":  cast(FloatArray, np.array([[-1, 0, 0], [ 0, 0, 1], [ 0, 1, 0]], dtype=float)),
    "FFS":  cast(FloatArray, np.array([[-1, 0, 0], [ 0, 0,-1], [ 0,-1, 0]], dtype=float)),
    "FFP":  cast(FloatArray, np.array([[ 1, 0, 0], [ 0, 0, 1], [ 0,-1, 0]], dtype=float)),
    "HFDL": cast(FloatArray, np.array([[ 0, 0,-1], [-1, 0, 0], [ 0, 1, 0]], dtype=float)),
    "HFDR": cast(FloatArray, np.array([[ 0, 0, 1], [ 1, 0, 0], [ 0, 1, 0]], dtype=float)),
    "FFDL": cast(FloatArray, np.array([[ 0, 0,-1], [ 1, 0, 0], [ 0,-1, 0]], dtype=float)),
    "FFDR": cast(FloatArray, np.array([[ 0, 0, 1], [-1, 0, 0], [ 0,-1, 0]], dtype=float)),
}


def _c_rot_z(row: FloatArray, alpha: float) -> None:
    """osh_vect_rot_z convention: x'=cos*x+sin*y, y'=-sin*x+cos*y (= Rz(-alpha))."""
    c, s = math.cos(alpha), math.sin(alpha)
    x, y = row[0], row[1]
    row[0] = c * x + s * y
    row[1] = -s * x + c * y


def _c_rot_y(row: FloatArray, alpha: float) -> None:
    """osh_vect_rot_y convention: x'=cos*x-sin*z, z'=sin*x+cos*z (= Ry(-alpha))."""
    c, s = math.cos(alpha), math.sin(alpha)
    x, z = row[0], row[2]
    row[0] = c * x - s * z
    row[2] = s * x + c * z


def parse_dcm_origin_from_geo(path: Path, dcm_path: Path | None) -> tuple[Path, FloatArray] | None:
    """Return (ct_dir, t_corner_cm) for the first matching DCM card in a geo.dat.

    DCM card format (current):
        DCM name ct_dir patient_pos gantry_deg couch_deg iso_x_mm iso_y_mm iso_z_mm
    where iso_{x,y,z}_mm is the treatment isocenter in DICOM LPS patient coords [mm].
    The CT corner in universe [cm] is computed using the same transform chain as the C parser.
    """
    pydicom = _import_pydicom()
    _f = r"([-+0-9.eE]+)"
    pattern = re.compile(
        r"^\s*DCM\s+\S+\s+(\S+)\s+(\S+)\s+" + _f + r"\s+" + _f + r"\s+" + _f + r"\s+" + _f + r"\s+" + _f
    )
    for line in path.read_text(encoding="utf-8").splitlines():
        match = pattern.match(line)
        if not match:
            continue
        geo_dcm_path = (path.parent / match.group(1)).resolve()
        if dcm_path is not None and geo_dcm_path != dcm_path.resolve():
            continue

        patient_pos = match.group(2).upper()
        gantry_rad = math.radians(float(match.group(3)))
        couch_rad = math.radians(float(match.group(4)))
        iso_mm = np.array([float(match.group(i)) for i in (5, 6, 7)], dtype=float)

        sample = next((p for p in sorted(geo_dcm_path.iterdir()) if p.is_file()), None)
        if sample is None:
            raise ValueError(f"{geo_dcm_path}: no DICOM files found")
        ds = pydicom.dcmread(str(sample), stop_before_pixels=True)
        spacing_cm = cast(FloatArray, np.array([
            0.1 * float(ds.PixelSpacing[1]),
            0.1 * float(ds.PixelSpacing[0]),
            0.1 * float(getattr(ds, "SliceThickness", 1.0)),
        ], dtype=float))
        ct_origin_mm = np.array([float(v) for v in ds.ImagePositionPatient], dtype=float)

        # Replicate C transform chain exactly (osh_geometry_parse.c):
        #   1. patient-position base rotation
        #   2. rot_z(-couch_rad, row)  — negated because rot_z is CW
        #   3. rot_y(gantry_rad, row)
        tb = _PP_BASE_ROTATION.get(patient_pos, _PP_BASE_ROTATION["HFS"]).copy()
        for row in tb:
            _c_rot_z(row, -couch_rad)
            _c_rot_y(row, gantry_rad)

        # iso_ct_local[j] = (iso_mm[j] - ct_origin_mm[j]) / 10 + spacing_cm[j] / 2
        iso_ct_local = (iso_mm - ct_origin_mm) * 0.1 + spacing_cm * 0.5

        # t_corner[k] = -sum_j tb[j][k] * iso_ct_local[j]
        t_corner = cast(FloatArray, -(tb.T @ iso_ct_local))
        return geo_dcm_path, t_corner
    return None


def _centers(lo: float, hi: float, n: int) -> FloatArray:
    edges = np.linspace(lo, hi, n + 1)
    return cast(FloatArray, 0.5 * (edges[:-1] + edges[1:]))


def _nearest_index(values: FloatArray, target: float) -> int:
    return int(np.argmin(np.abs(values - target)))


def extract_planes(
    ct: CtVolume,
    bdo: BdoData,
    page_idx: int,
    plane: str = "auto",
    cut_pos: float | None = None,
) -> tuple[FloatArray, FloatArray, tuple[float, float, float, float], str]:
    """Return (ct_plane, score_plane, ct_extent, plane_name).

    plane: 'XZ', 'XY', 'YZ', or 'auto' (picks the singleton axis if present,
           otherwise defaults to XZ through the centre of the scoring grid).
    cut_pos: position along the cut axis in cm; defaults to grid centre.
    """
    nx, ny, nz = (int(v) for v in bdo.n)
    score = bdo.pages[page_idx].reshape((nz, ny, nx))

    x_ct = cast(FloatArray, ct.origin_cm[0] + ct.spacing_cm[0] * (np.arange(ct.data_hu.shape[2]) + 0.5))
    y_ct = cast(FloatArray, ct.origin_cm[1] + ct.spacing_cm[1] * (np.arange(ct.data_hu.shape[1]) + 0.5))
    z_ct = cast(FloatArray, ct.origin_cm[2] + ct.spacing_cm[2] * (np.arange(ct.data_hu.shape[0]) + 0.5))

    # Resolve 'auto': prefer singleton axis, else default to XZ.
    if plane == "auto":
        singleton_axes = [i for i, v in enumerate((nx, ny, nz)) if v == 1]
        if len(singleton_axes) == 1:
            plane = ["YZ", "XZ", "XY"][singleton_axes[0]]
        else:
            plane = "XZ"

    x_lo, x_hi = float(bdo.p[0]), float(bdo.q[0])
    y_lo, y_hi = float(bdo.p[1]), float(bdo.q[1])
    z_lo, z_hi = float(bdo.p[2]), float(bdo.q[2])

    if plane == "XZ":
        y_centers = _centers(y_lo, y_hi, ny)
        y0 = cut_pos if cut_pos is not None else 0.5 * (y_lo + y_hi)
        iy_score = _nearest_index(y_centers, y0)
        iy_ct = _nearest_index(y_ct, y_centers[iy_score])
        ct_plane = cast(FloatArray, ct.data_hu[:, iy_ct, :])
        score_plane = cast(FloatArray, score[:, iy_score, :])
        return ct_plane, score_plane, (x_ct[0], x_ct[-1], z_ct[0], z_ct[-1]), "XZ"

    if plane == "YZ":
        x_centers = _centers(x_lo, x_hi, nx)
        x0 = cut_pos if cut_pos is not None else 0.5 * (x_lo + x_hi)
        ix_score = _nearest_index(x_centers, x0)
        ix_ct = _nearest_index(x_ct, x_centers[ix_score])
        ct_plane = cast(FloatArray, ct.data_hu[:, :, ix_ct])
        score_plane = cast(FloatArray, score[:, :, ix_score])
        return ct_plane, score_plane, (y_ct[0], y_ct[-1], z_ct[0], z_ct[-1]), "YZ"

    # XY
    z_centers = _centers(z_lo, z_hi, nz)
    z0 = cut_pos if cut_pos is not None else 0.5 * (z_lo + z_hi)
    iz_score = _nearest_index(z_centers, z0)
    iz_ct = _nearest_index(z_ct, z_centers[iz_score])
    ct_plane = cast(FloatArray, ct.data_hu[iz_ct, :, :])
    score_plane = cast(FloatArray, score[iz_score, :, :])
    return ct_plane, score_plane, (x_ct[0], x_ct[-1], y_ct[0], y_ct[-1]), "XY"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("dicom", type=Path,
                        help="DICOM CT directory or single CT file")
    parser.add_argument("scoring", type=Path,
                        help="scoring output: BDO2019 file or DICOM RTDOSE (.dcm)")
    parser.add_argument("-o", "--output", type=Path,
                        help="write plot to this image file instead of showing interactively")
    parser.add_argument("--geo", type=Path,
                        help="geo.dat; DCM card is used to infer the CT simulation origin")
    parser.add_argument("--ct-origin", nargs=3, type=float, metavar=("X", "Y", "Z"),
                        help="CT voxel-corner origin in cm (overrides geo.dat and DICOM metadata)")
    parser.add_argument("--page", type=int, default=0,
                        help="scoring page index to plot (default: 0)")
    parser.add_argument("--plane", choices=["XZ", "XY", "YZ", "auto"], default="auto",
                        help="projection plane (default: auto — singleton axis, else XZ)")
    parser.add_argument("--cut", type=float, default=None, metavar="CM",
                        help="position along the cut axis in cm (default: grid centre)")
    parser.add_argument("--window", nargs=2, type=float, default=(-1000.0, 1000.0), metavar=("LO", "HI"),
                        help="CT HU window (default: -1000 1000)")
    parser.add_argument("--alpha", type=float, default=0.6,
                        help="score overlay opacity (default: 0.6)")
    parser.add_argument("--log", action="store_true",
                        help="plot positive score values on a log10 scale")
    args = parser.parse_args()

    dicom_path = cast(Path, args.dicom)
    origin = cast(FloatArray, np.asarray(args.ct_origin, dtype=np.float64)) if args.ct_origin is not None else None
    if args.geo is not None and origin is None:
        parsed = parse_dcm_origin_from_geo(args.geo, dicom_path)
        if parsed is not None:
            dicom_path, origin = parsed

    bdo = read_scoring(args.scoring)
    ct = read_ct_series(dicom_path, origin)
    if args.page < 0 or args.page >= len(bdo.pages):
        raise ValueError(f"{args.scoring}: page {args.page} out of range, file has {len(bdo.pages)} page(s)")

    ct_plane, score_plane, ct_extent, plane = extract_planes(ct, bdo, args.page, args.plane, args.cut)

    if plane == "XZ":
        score_extent = (float(bdo.p[0]), float(bdo.q[0]), float(bdo.p[2]), float(bdo.q[2]))
    elif plane == "YZ":
        score_extent = (float(bdo.p[1]), float(bdo.q[1]), float(bdo.p[2]), float(bdo.q[2]))
    else:
        score_extent = (float(bdo.p[0]), float(bdo.q[0]), float(bdo.p[1]), float(bdo.q[1]))

    overlay = np.ma.masked_where(score_plane <= 0.0, score_plane)
    if args.log:
        overlay = np.ma.log10(overlay)

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.imshow(ct_plane, cmap="gray", origin="lower", extent=ct_extent,
              vmin=args.window[0], vmax=args.window[1], aspect="equal")
    im = ax.imshow(overlay, cmap="inferno", origin="lower", extent=score_extent,
                   alpha=args.alpha, aspect="equal")
    units = bdo.page_units[args.page] if args.page < len(bdo.page_units) else "arb"
    label = f"log10(score / {units})" if args.log else f"score ({units})"
    fig.colorbar(im, ax=ax, label=label)
    ax.set_title(f"{bdo.geo_name or args.scoring.name} page {args.page} on CT ({plane})")
    ax.set_xlabel(f"{plane[0]} [cm]")
    ax.set_ylabel(f"{plane[1]} [cm]")
    ax.set_xlim(score_extent[0], score_extent[1])
    ax.set_ylim(score_extent[2], score_extent[3])
    fig.tight_layout()

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.output, dpi=160)
        print(f"saved: {args.output}")
    else:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
