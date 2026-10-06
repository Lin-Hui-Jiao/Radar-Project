#!/usr/bin/env python3
"""Run BVH-RT occlusion checks from a fixed-height candidate CSV.

Input candidate CSV must contain local_x, local_y, local_z. The transmitter
is supplied as lon/lat/alt and converted to the same local coordinate system
used by the generated candidate file.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from bvh_rt_multiobj_experiment import (  # noqa: E402
    build_scene,
    cast_visibility_batch,
    current_rss_gb,
    discover_obj_paths,
    get_open3d,
    scan_obj_bbox,
    write_csv,
)


DEFAULT_RADAR_LON = 114.1670
DEFAULT_RADAR_LAT = 22.2806
DEFAULT_RADAR_ALT = 164.0
DEFAULT_INDEX_RANGE_X = -1975.0
DEFAULT_INDEX_RANGE_Y = 43.0
DEFAULT_MIN_X = 835000.0
DEFAULT_MIN_Y = 815500.0


def import_pyproj():
    try:
        from pyproj import Transformer
    except ImportError as exc:
        raise SystemExit("pyproj is required: python -m pip install pyproj") from exc
    return Transformer


def now_s() -> float:
    return time.perf_counter()


def local_from_lonlat(
    lon: float,
    lat: float,
    alt: float,
    source_crs: str,
    projected_crs: str,
    min_x: float,
    min_y: float,
    index_range_x: float,
    index_range_y: float,
) -> np.ndarray:
    Transformer = import_pyproj()
    transformer = Transformer.from_crs(source_crs, projected_crs, always_xy=True)
    x2326, y2326 = transformer.transform(lon, lat)
    return np.array(
        [
            index_range_x + (x2326 - min_x),
            index_range_y + (y2326 - min_y),
            alt,
        ],
        dtype=np.float64,
    )


def load_candidate_csv(path: Path, chunk_rows: int) -> tuple[np.ndarray, list[str]]:
    if chunk_rows <= 0:
        raise ValueError("--csv-chunk-rows must be positive")

    chunks: list[np.ndarray] = []
    fieldnames: list[str] = []
    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"Candidate CSV has no header: {path}")
        fieldnames = list(reader.fieldnames)
        required = {"local_x", "local_y", "local_z"}
        missing = required.difference(fieldnames)
        if missing:
            raise ValueError(f"Candidate CSV is missing columns: {sorted(missing)}")

        buf = np.empty((chunk_rows, 3), dtype=np.float64)
        count = 0
        for row in reader:
            try:
                buf[count, 0] = float(row["local_x"])
                buf[count, 1] = float(row["local_y"])
                buf[count, 2] = float(row["local_z"])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Invalid local coordinate row near index {sum(c.shape[0] for c in chunks) + count}") from exc
            count += 1
            if count == chunk_rows:
                chunks.append(buf.copy())
                count = 0
        if count:
            chunks.append(buf[:count].copy())

    if not chunks:
        raise ValueError(f"Candidate CSV contains no points: {path}")
    points = np.vstack(chunks)
    return np.ascontiguousarray(points, dtype=np.float64), fieldnames


def write_result_csv(
    input_csv: Path,
    output_csv: Path,
    visible: np.ndarray,
    blocked_column: str,
    visible_column: str,
) -> None:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with input_csv.open("r", newline="", encoding="utf-8") as f_in, output_csv.open(
        "w", newline="", encoding="utf-8"
    ) as f_out:
        reader = csv.DictReader(f_in)
        if reader.fieldnames is None:
            raise ValueError(f"Candidate CSV has no header: {input_csv}")
        fieldnames = list(reader.fieldnames)
        for col in [visible_column, blocked_column]:
            if col not in fieldnames:
                fieldnames.append(col)
        writer = csv.DictWriter(f_out, fieldnames=fieldnames)
        writer.writeheader()
        for idx, row in enumerate(reader):
            row[visible_column] = "1" if bool(visible[idx]) else "0"
            row[blocked_column] = "0" if bool(visible[idx]) else "1"
            writer.writerow(row)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Read fixed-height candidate CSV and run BVH-RT OBJ occlusion checks."
    )
    parser.add_argument("--candidates", required=True, help="Candidate CSV with local_x,local_y,local_z columns.")
    parser.add_argument("--obj", action="append", default=[], help="OBJ file, directory, or glob. May be repeated.")
    parser.add_argument("--obj-dir", help="Directory containing OBJ tiles.")
    parser.add_argument("--obj-glob", default="*.obj")
    parser.add_argument("--radar-lon", type=float, default=DEFAULT_RADAR_LON)
    parser.add_argument("--radar-lat", type=float, default=DEFAULT_RADAR_LAT)
    parser.add_argument("--radar-alt", type=float, default=DEFAULT_RADAR_ALT)
    parser.add_argument("--source-crs", default="EPSG:4326")
    parser.add_argument("--projected-crs", default="EPSG:2326")
    parser.add_argument("--min-x", type=float, default=DEFAULT_MIN_X)
    parser.add_argument("--min-y", type=float, default=DEFAULT_MIN_Y)
    parser.add_argument("--index-range-x", type=float, default=DEFAULT_INDEX_RANGE_X)
    parser.add_argument("--index-range-y", type=float, default=DEFAULT_INDEX_RANGE_Y)
    parser.add_argument("--eps", type=float, default=0.01)
    parser.add_argument("--ray-chunk-size", type=int, default=500_000)
    parser.add_argument("--csv-chunk-rows", type=int, default=200_000)
    parser.add_argument("--build-threads", type=int, default=0)
    parser.add_argument("--query-threads", type=int, default=0)
    parser.add_argument("--mesh-offset", nargs=3, type=float, help="Manual offset subtracted from mesh/points.")
    parser.add_argument("--zero-mesh-offset", action="store_true", help="Do not apply Open3D numeric-stability offset.")
    parser.add_argument("--output-csv", default="out/bvh_rt_candidate_results.csv")
    parser.add_argument("--labels-output", default="out/bvh_rt_candidate_visible.npy")
    parser.add_argument("--summary-output", default="out/bvh_rt_candidate_summary.csv")
    parser.add_argument("--mesh-report", default="out/bvh_rt_candidate_mesh_report.csv")
    parser.add_argument("--visible-column", default="bvh_visible")
    parser.add_argument("--blocked-column", default="bvh_blocked")
    parser.add_argument("--progress-every", type=int, default=0)
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    candidate_csv = Path(args.candidates)

    load_start = now_s()
    targets, _ = load_candidate_csv(candidate_csv, args.csv_chunk_rows)
    candidate_load_time_s = now_s() - load_start
    print(f"Loaded {targets.shape[0]} candidate points from {candidate_csv}")

    radar_local = local_from_lonlat(
        lon=args.radar_lon,
        lat=args.radar_lat,
        alt=args.radar_alt,
        source_crs=args.source_crs,
        projected_crs=args.projected_crs,
        min_x=args.min_x,
        min_y=args.min_y,
        index_range_x=args.index_range_x,
        index_range_y=args.index_range_y,
    )
    print(f"Radar local xyz: {radar_local.tolist()}")

    o3d = get_open3d()
    obj_paths = discover_obj_paths(args.obj, args.obj_dir, args.obj_glob)
    print(f"Discovered {len(obj_paths)} OBJ files")

    bbox_scan_time_s = 0.0
    bbox_vertex_count = 0
    bbox_min = np.array([float("nan")] * 3, dtype=np.float64)
    bbox_max = np.array([float("nan")] * 3, dtype=np.float64)
    if args.mesh_offset is not None:
        mesh_offset = np.array(args.mesh_offset, dtype=np.float64)
    elif args.zero_mesh_offset:
        mesh_offset = np.zeros(3, dtype=np.float64)
    else:
        bbox_min, bbox_max, bbox_vertex_count, bbox_scan_time_s = scan_obj_bbox(obj_paths)
        mesh_offset = (bbox_min + bbox_max) * 0.5
    print(f"Open3D numeric offset: {mesh_offset.tolist()}")

    rss_before_build_gb = current_rss_gb()
    scene, build_stats, mesh_rows = build_scene(o3d, obj_paths, mesh_offset, args.build_threads)
    rss_after_build_gb = current_rss_gb()
    write_csv(Path(args.mesh_report), mesh_rows)

    compute_start = now_s()
    visible, query_stats = cast_visibility_batch(
        scene=scene,
        tx=radar_local,
        targets=targets,
        offset=mesh_offset,
        eps=args.eps,
        chunk_size=args.ray_chunk_size,
        query_threads=args.query_threads,
        progress_every=args.progress_every,
    )
    measured_compute_time_s = now_s() - compute_start
    core_compute_time_s = float(query_stats.get("query_time_s", measured_compute_time_s))
    cast_rays_time_s = float(query_stats.get("cast_rays_time_s", float("nan")))
    ray_prepare_and_classify_time_s = float(
        query_stats.get(
            "ray_prepare_and_classify_time_s",
            max(0.0, core_compute_time_s - cast_rays_time_s)
            if math.isfinite(cast_rays_time_s)
            else float("nan"),
        )
    )
    cast_rays_throughput_qps = float(
        query_stats.get(
            "cast_rays_throughput_qps",
            float(query_stats["num_queries"]) / cast_rays_time_s
            if cast_rays_time_s > 0.0
            else float("nan"),
        )
    )
    throughput_qps = float(
        query_stats.get(
            "throughput_qps",
            float(query_stats.get("num_queries", len(targets))) / core_compute_time_s
            if core_compute_time_s > 0.0
            else float("nan"),
        )
    )

    labels_path = Path(args.labels_output)
    labels_path.parent.mkdir(parents=True, exist_ok=True)
    np.save(labels_path, visible)
    write_result_csv(
        input_csv=candidate_csv,
        output_csv=Path(args.output_csv),
        visible=visible,
        visible_column=args.visible_column,
        blocked_column=args.blocked_column,
    )

    row: dict[str, object] = {
        "method": "BVH-RT-Open3D",
        "candidate_file": str(candidate_csv.resolve(strict=False)),
        "output_csv": str(Path(args.output_csv).resolve(strict=False)),
        "labels_output": str(labels_path.resolve(strict=False)),
        "num_obj_files": build_stats["num_obj_files"],
        "total_vertices": build_stats["total_vertices"],
        "total_triangles": build_stats["total_triangles"],
        "radar_lon": args.radar_lon,
        "radar_lat": args.radar_lat,
        "radar_alt": args.radar_alt,
        "radar_local_x": radar_local[0],
        "radar_local_y": radar_local[1],
        "radar_local_z": radar_local[2],
        "min_x": args.min_x,
        "min_y": args.min_y,
        "index_range_x": args.index_range_x,
        "index_range_y": args.index_range_y,
        "bbox_min_x": bbox_min[0],
        "bbox_min_y": bbox_min[1],
        "bbox_min_z": bbox_min[2],
        "bbox_max_x": bbox_max[0],
        "bbox_max_y": bbox_max[1],
        "bbox_max_z": bbox_max[2],
        "bbox_vertex_count": bbox_vertex_count,
        "mesh_offset_x": mesh_offset[0],
        "mesh_offset_y": mesh_offset[1],
        "mesh_offset_z": mesh_offset[2],
        "eps": args.eps,
        "ray_chunk_size": args.ray_chunk_size,
        "build_threads": args.build_threads,
        "query_threads": args.query_threads,
        "candidate_load_time_s": candidate_load_time_s,
        "bbox_scan_time_s": bbox_scan_time_s,
        "mesh_read_time_s": build_stats["mesh_read_time_s"],
        "mesh_clean_shift_time_s": build_stats["mesh_clean_shift_time_s"],
        "scene_add_time_s": build_stats["scene_add_time_s"],
        "scene_build_time_s": build_stats["scene_build_time_s"],
        "rss_before_build_gb": rss_before_build_gb,
        "rss_after_build_gb": rss_after_build_gb,
        "num_queries": query_stats["num_queries"],
        "core_compute_time_s": core_compute_time_s,
        "query_time_s": core_compute_time_s,
        "cast_rays_time_s": cast_rays_time_s,
        "ray_prepare_and_classify_time_s": ray_prepare_and_classify_time_s,
        "throughput_qps": throughput_qps,
        "cast_rays_throughput_qps": cast_rays_throughput_qps,
        "num_visible": query_stats["num_visible"],
        "num_blocked": query_stats["num_blocked"],
        "num_invalid_targets": query_stats["num_invalid_targets"],
        "num_near_tx_targets": query_stats["num_near_tx_targets"],
        "visible_ratio": query_stats["visible_ratio"],
        "peak_rss_gb": query_stats["peak_rss_gb"],
    }
    write_csv(Path(args.summary_output), [row])
    meta_path = Path(args.summary_output).with_suffix(Path(args.summary_output).suffix + ".json")
    with meta_path.open("w", encoding="utf-8") as f:
        json.dump(row, f, ensure_ascii=False, indent=2)

    print(f"Result CSV written to {args.output_csv}")
    print(f"Labels written to {labels_path}")
    print(f"Summary written to {args.summary_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
