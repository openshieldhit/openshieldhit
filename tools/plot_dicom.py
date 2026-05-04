#!/usr/bin/env python3
"""Overlay an OpenShieldHIT BDO scoring plane on a DICOM CT slice."""

from __future__ import annotations

import argparse
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
    p: FloatArray
    q: FloatArray
    n: IntArray
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


def parse_dcm_origin_from_geo(path: Path, dcm_path: Path | None) -> tuple[Path, FloatArray] | None:
    pydicom = _import_pydicom()
    pattern = re.compile(r"^\s*DCM\s+\S+\s+(\S+)\s+\S+\s+\S+\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)")
    for line in path.read_text(encoding="utf-8").splitlines():
        match = pattern.match(line)
        if not match:
            continue
        geo_dcm_path = (path.parent / match.group(1)).resolve()
        if dcm_path is not None and geo_dcm_path != dcm_path.resolve():
            continue
        tx, ty, tz = (float(match.group(i)) for i in (2, 3, 4))
        sample = next((p for p in sorted(geo_dcm_path.iterdir()) if p.is_file()), None)
        if sample is None:
            raise ValueError(f"{geo_dcm_path}: no DICOM files found")
        ds = pydicom.dcmread(str(sample), stop_before_pixels=True)
        dx = 0.1 * float(ds.PixelSpacing[1])
        dy = 0.1 * float(ds.PixelSpacing[0])
        dz = 0.1 * float(getattr(ds, "SliceThickness", 1.0))
        origin_cm = cast(FloatArray, np.asarray([tx - 0.5 * dx, ty - 0.5 * dy, tz - 0.5 * dz], dtype=np.float64))
        return geo_dcm_path, origin_cm
    return None


def _centers(lo: float, hi: float, n: int) -> FloatArray:
    edges = np.linspace(lo, hi, n + 1)
    return cast(FloatArray, 0.5 * (edges[:-1] + edges[1:]))


def _nearest_index(values: FloatArray, target: float) -> int:
    return int(np.argmin(np.abs(values - target)))


def extract_planes(ct: CtVolume, bdo: BdoData, page_idx: int) -> tuple[FloatArray, FloatArray, tuple[float, float, float, float], str]:
    nx, ny, nz = (int(v) for v in bdo.n)
    score = bdo.pages[page_idx].reshape((nz, ny, nx))
    singleton_axes = [i for i, n in enumerate((nx, ny, nz)) if n == 1]
    if len(singleton_axes) != 1:
        raise ValueError(f"{bdo.path}: expected one singleton scoring axis for a 2D overlay, got n={bdo.n.tolist()}")

    x_edges = (float(bdo.p[0]), float(bdo.q[0]))
    y_edges = (float(bdo.p[1]), float(bdo.q[1]))
    z_edges = (float(bdo.p[2]), float(bdo.q[2]))
    x_ct = cast(FloatArray, ct.origin_cm[0] + ct.spacing_cm[0] * (np.arange(ct.data_hu.shape[2]) + 0.5))
    y_ct = cast(FloatArray, ct.origin_cm[1] + ct.spacing_cm[1] * (np.arange(ct.data_hu.shape[1]) + 0.5))
    z_ct = cast(FloatArray, ct.origin_cm[2] + ct.spacing_cm[2] * (np.arange(ct.data_hu.shape[0]) + 0.5))

    if singleton_axes[0] == 1:
        y0 = _centers(*y_edges, ny)[0]
        ct_plane = cast(FloatArray, ct.data_hu[:, _nearest_index(y_ct, y0), :])
        score_plane = cast(FloatArray, score[:, 0, :])
        return ct_plane, score_plane, (x_ct[0], x_ct[-1], z_ct[0], z_ct[-1]), "XZ"
    if singleton_axes[0] == 0:
        x0 = _centers(*x_edges, nx)[0]
        ct_plane = cast(FloatArray, ct.data_hu[:, :, _nearest_index(x_ct, x0)])
        score_plane = cast(FloatArray, score[:, :, 0])
        return ct_plane, score_plane, (y_ct[0], y_ct[-1], z_ct[0], z_ct[-1]), "YZ"

    z0 = _centers(*z_edges, nz)[0]
    ct_plane = cast(FloatArray, ct.data_hu[_nearest_index(z_ct, z0), :, :])
    score_plane = cast(FloatArray, score[0, :, :])
    return ct_plane, score_plane, (x_ct[0], x_ct[-1], y_ct[0], y_ct[-1]), "XY"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dicom", type=Path, help="DICOM CT file or series directory")
    parser.add_argument("bdo", type=Path, help="OpenShieldHIT BDO2019 scoring output")
    parser.add_argument("-o", "--output", type=Path, help="write plot to this image instead of showing it")
    parser.add_argument("--geo", type=Path, help="optional geo.dat; DCM card is used to infer CT simulation origin")
    parser.add_argument("--ct-origin", nargs=3, type=float, metavar=("X", "Y", "Z"), help="CT voxel-corner origin in cm")
    parser.add_argument("--page", type=int, default=0, help="BDO page index to plot")
    parser.add_argument("--window", nargs=2, type=float, default=(-1000.0, 1000.0), metavar=("LO", "HI"))
    parser.add_argument("--alpha", type=float, default=0.6, help="score overlay opacity")
    parser.add_argument("--log", action="store_true", help="plot positive score values on a log10 scale")
    args = parser.parse_args()

    dicom_path = cast(Path, args.dicom)
    origin = cast(FloatArray, np.asarray(args.ct_origin, dtype=np.float64)) if args.ct_origin is not None else None
    if args.geo is not None and origin is None:
        parsed = parse_dcm_origin_from_geo(args.geo, dicom_path)
        if parsed is not None:
            dicom_path, origin = parsed

    bdo = read_bdo2019(args.bdo)
    ct = read_ct_series(dicom_path, origin)
    if args.page < 0 or args.page >= len(bdo.pages):
        raise ValueError(f"{args.bdo}: page {args.page} out of range, file has {len(bdo.pages)} page(s)")

    ct_plane, score_plane, ct_extent, plane = extract_planes(ct, bdo, args.page)
    score_extent = (float(bdo.p[0]), float(bdo.q[0]), float(bdo.p[2]), float(bdo.q[2]))
    if plane == "YZ":
        score_extent = (float(bdo.p[1]), float(bdo.q[1]), float(bdo.p[2]), float(bdo.q[2]))
    elif plane == "XY":
        score_extent = (float(bdo.p[0]), float(bdo.q[0]), float(bdo.p[1]), float(bdo.q[1]))

    overlay = np.ma.masked_where(score_plane <= 0.0, score_plane)
    if args.log:
        overlay = np.ma.log10(overlay)

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.imshow(ct_plane, cmap="gray", origin="lower", extent=ct_extent, vmin=args.window[0], vmax=args.window[1], aspect="equal")
    im = ax.imshow(overlay, cmap="inferno", origin="lower", extent=score_extent, alpha=args.alpha, aspect="equal")
    units = bdo.page_units[args.page] if args.page < len(bdo.page_units) else "arb"
    label = f"log10(score / {units})" if args.log else f"score ({units})"
    fig.colorbar(im, ax=ax, label=label)
    ax.set_title(f"{bdo.geo_name or args.bdo.name} page {args.page} on CT ({plane})")
    ax.set_xlabel(f"{plane[0]} [cm]")
    ax.set_ylabel(f"{plane[1]} [cm]")
    ax.set_xlim(score_extent[0], score_extent[1])
    ax.set_ylim(score_extent[2], score_extent[3])
    fig.tight_layout()

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.output, dpi=160)
    else:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
