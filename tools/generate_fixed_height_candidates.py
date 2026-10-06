#!/usr/bin/env python3
"""Generate fixed-height local-coordinate candidate points from a lon/lat dataset range.

The output CSV is intended to be the common input for both VRPF and BVH-RT.
Rows are sampled on the local projected coordinate grid, filtered so their
inverse-projected lon/lat stays inside the requested dataset extent, and then
written mainly as local x/y/z coordinates.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from pathlib import Path


DEFAULT_MIN_LON = 114.164571
DEFAULT_MAX_LON = 114.169422
DEFAULT_MIN_LAT = 22.278365
DEFAULT_MAX_LAT = 22.282880

DEFAULT_INDEX_RANGE_X = -1975.0
DEFAULT_INDEX_RANGE_Y = 43.0
DEFAULT_MIN_X = 835000.0
DEFAULT_MIN_Y = 815500.0
DEFAULT_RADAR_LON = 114.1670
DEFAULT_RADAR_LAT = 22.2806
DEFAULT_RADAR_ALT = 80.0
DEFAULT_MIN_ALT = 0.0
DEFAULT_MAX_ALT = 100.0


def import_pyproj():
    try:
        from pyproj import Transformer
    except ImportError as exc:
        raise SystemExit(
            "pyproj is required. Install it in the Python environment used for "
            "this script, for example: python -m pip install pyproj"
        ) from exc
    return Transformer


def frange(start: float, stop: float, step: float):
    value = start
    # Include the far edge when it lands exactly on the grid within tolerance.
    while value <= stop + step * 1e-9:
        yield value
        value += step


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a fixed-height candidate CSV in the project local coordinate system."
    )
    parser.add_argument("--min-lon", type=float, default=DEFAULT_MIN_LON)
    parser.add_argument("--max-lon", type=float, default=DEFAULT_MAX_LON)
    parser.add_argument("--min-lat", type=float, default=DEFAULT_MIN_LAT)
    parser.add_argument("--max-lat", type=float, default=DEFAULT_MAX_LAT)
    parser.add_argument("--radar-lon", type=float, default=DEFAULT_RADAR_LON)
    parser.add_argument("--radar-lat", type=float, default=DEFAULT_RADAR_LAT)
    parser.add_argument("--radar-alt", type=float, default=DEFAULT_RADAR_ALT)
    parser.add_argument("--min-alt", type=float, default=DEFAULT_MIN_ALT)
    parser.add_argument("--max-alt", type=float, default=DEFAULT_MAX_ALT)
    parser.add_argument("--fixed-height", "--z", dest="fixed_height", type=float, required=True)
    parser.add_argument("--spacing", type=float, default=1.0, help="Sampling spacing in local metres.")
    parser.add_argument("--output", required=True, help="Output CSV path.")
    parser.add_argument("--meta-output", help="Optional metadata JSON path. Default: <output>.meta.json")
    parser.add_argument("--source-crs", default="EPSG:4326")
    parser.add_argument("--projected-crs", default="EPSG:2326")
    parser.add_argument("--min-x", type=float, default=DEFAULT_MIN_X, help="EPSG:2326 x offset used by the dataset.")
    parser.add_argument("--min-y", type=float, default=DEFAULT_MIN_Y, help="EPSG:2326 y offset used by the dataset.")
    parser.add_argument("--index-range-x", type=float, default=DEFAULT_INDEX_RANGE_X)
    parser.add_argument("--index-range-y", type=float, default=DEFAULT_INDEX_RANGE_Y)
    parser.add_argument("--chunk-rows", type=int, default=200_000)
    parser.add_argument(
        "--local-only",
        action="store_true",
        help="Write only point_id,local_x,local_y,local_z. By default lon/lat are also included.",
    )
    parser.add_argument("--precision", type=int, default=6, help="Decimal digits for local coordinates.")
    parser.add_argument("--lonlat-precision", type=int, default=10)
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
    if args.min_lon >= args.max_lon:
        raise ValueError("--min-lon must be less than --max-lon")
    if args.min_lat >= args.max_lat:
        raise ValueError("--min-lat must be less than --max-lat")
    if args.spacing <= 0:
        raise ValueError("--spacing must be positive")
    if args.min_alt > args.max_alt:
        raise ValueError("--min-alt must be less than or equal to --max-alt")
    if not (args.min_alt <= args.fixed_height <= args.max_alt):
        raise ValueError("--fixed-height must lie within [--min-alt, --max-alt]")
    if args.chunk_rows <= 0:
        raise ValueError("--chunk-rows must be positive")
    if not math.isfinite(args.fixed_height):
        raise ValueError("--fixed-height must be finite")
    if not all(math.isfinite(v) for v in [args.radar_lon, args.radar_lat, args.radar_alt]):
        raise ValueError("Radar parameters must be finite")


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    validate_args(args)

    Transformer = import_pyproj()
    to_projected = Transformer.from_crs(args.source_crs, args.projected_crs, always_xy=True)
    to_lonlat = Transformer.from_crs(args.projected_crs, args.source_crs, always_xy=True)

    corners_lon = [args.min_lon, args.max_lon, args.max_lon, args.min_lon]
    corners_lat = [args.min_lat, args.min_lat, args.max_lat, args.max_lat]
    corners_x2326, corners_y2326 = to_projected.transform(corners_lon, corners_lat)

    local_x_corners = [args.index_range_x + (x - args.min_x) for x in corners_x2326]
    local_y_corners = [args.index_range_y + (y - args.min_y) for y in corners_y2326]
    local_min_x = math.floor(min(local_x_corners) / args.spacing) * args.spacing
    local_max_x = math.ceil(max(local_x_corners) / args.spacing) * args.spacing
    local_min_y = math.floor(min(local_y_corners) / args.spacing) * args.spacing
    local_max_y = math.ceil(max(local_y_corners) / args.spacing) * args.spacing

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    start_time = time.perf_counter()
    total_candidates = 0
    kept = 0

    if args.local_only:
        fieldnames = ["point_id", "local_x", "local_y", "local_z"]
    else:
        fieldnames = ["point_id", "local_x", "local_y", "local_z", "lon", "lat"]

    local_fmt = f"%.{args.precision}f"
    lonlat_fmt = f"%.{args.lonlat_precision}f"

    with output.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        buffer: list[dict[str, str | int]] = []

        for local_y in frange(local_min_y, local_max_y, args.spacing):
            for local_x in frange(local_min_x, local_max_x, args.spacing):
                total_candidates += 1
                x2326 = local_x - args.index_range_x + args.min_x
                y2326 = local_y - args.index_range_y + args.min_y
                lon, lat = to_lonlat.transform(x2326, y2326)
                if not (args.min_lon <= lon <= args.max_lon and args.min_lat <= lat <= args.max_lat):
                    continue

                row: dict[str, str | int] = {
                    "point_id": kept,
                    "local_x": local_fmt % local_x,
                    "local_y": local_fmt % local_y,
                    "local_z": local_fmt % args.fixed_height,
                }
                if not args.local_only:
                    row["lon"] = lonlat_fmt % lon
                    row["lat"] = lonlat_fmt % lat
                buffer.append(row)
                kept += 1

                if len(buffer) >= args.chunk_rows:
                    writer.writerows(buffer)
                    buffer.clear()

        if buffer:
            writer.writerows(buffer)

    elapsed = time.perf_counter() - start_time
    meta = {
        "candidate_csv": str(output),
        "coordinate_system": {
            "source_crs": args.source_crs,
            "projected_crs": args.projected_crs,
            "local_x_formula": "index_range_x + (projected_x - min_x)",
            "local_y_formula": "index_range_y + (projected_y - min_y)",
            "z_unit": "metre",
        },
        "dataset": {
            "min_lon": args.min_lon,
            "max_lon": args.max_lon,
            "min_lat": args.min_lat,
            "max_lat": args.max_lat,
            "min_alt": args.min_alt,
            "max_alt": args.max_alt,
            "min_x": args.min_x,
            "min_y": args.min_y,
            "index_range_x": args.index_range_x,
            "index_range_y": args.index_range_y,
        },
        "radar": {
            "lon": args.radar_lon,
            "lat": args.radar_lat,
            "alt": args.radar_alt,
        },
        "slice": {
            "fixed_height": args.fixed_height,
            "spacing": args.spacing,
            "local_min_x": local_min_x,
            "local_max_x": local_max_x,
            "local_min_y": local_min_y,
            "local_max_y": local_max_y,
            "kept_points": kept,
            "candidate_grid_points_before_lonlat_filter": total_candidates,
        },
    }
    meta_output = Path(args.meta_output) if args.meta_output else output.with_suffix(output.suffix + ".meta.json")
    with meta_output.open("w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)

    print("Generated fixed-height candidate CSV")
    print(f"  output: {output}")
    print(f"  metadata: {meta_output}")
    print(f"  lon range: [{args.min_lon}, {args.max_lon}]")
    print(f"  lat range: [{args.min_lat}, {args.max_lat}]")
    print(f"  radar lon/lat/alt: {args.radar_lon}, {args.radar_lat}, {args.radar_alt}")
    print(f"  local bbox sampled: x=[{local_min_x}, {local_max_x}], y=[{local_min_y}, {local_max_y}]")
    print(f"  spacing: {args.spacing} m")
    print(f"  fixed height: {args.fixed_height}")
    print(f"  kept points: {kept}")
    print(f"  candidate grid points before lon/lat filter: {total_candidates}")
    print(f"  elapsed: {elapsed:.3f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
