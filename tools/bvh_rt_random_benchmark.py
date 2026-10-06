#!/usr/bin/env python3
"""Benchmark BVH-RT visibility queries on randomly sampled points.

The benchmark samples points inside a lon/lat region in batches, converts
them to the local coordinate system, and times only the BVH visibility query.
Random sampling, coordinate projection, OBJ loading, and scene construction are
reported separately and are not included in core_compute_time_s.
"""

from __future__ import annotations

import argparse
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


DEFAULT_MIN_LON = 114.164571
DEFAULT_MAX_LON = 114.169422
DEFAULT_MIN_LAT = 22.278365
DEFAULT_MAX_LAT = 22.282880
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


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be a positive integer")
    return parsed


def nonnegative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("value must be non-negative")
    return parsed


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
    x_projected, y_projected = transformer.transform(lon, lat)
    return np.array(
        [
            index_range_x + (x_projected - min_x),
            index_range_y + (y_projected - min_y),
            alt,
        ],
        dtype=np.float64,
    )


def generate_random_target_batch(
    rng: np.random.Generator,
    batch_size: int,
    transformer,
    args: argparse.Namespace,
) -> tuple[np.ndarray, dict[str, float]]:
    random_start = now_s()
    lons = rng.uniform(args.min_lon, args.max_lon, batch_size)
    lats = rng.uniform(args.min_lat, args.max_lat, batch_size)
    if args.random_alt:
        zs = rng.uniform(args.min_alt, args.max_alt, batch_size)
    else:
        zs = np.full(batch_size, args.fixed_height, dtype=np.float64)
    random_generation_time_s = now_s() - random_start

    projection_start = now_s()
    projected_x, projected_y = transformer.transform(lons, lats)
    local_x = args.index_range_x + (np.asarray(projected_x, dtype=np.float64) - args.min_x)
    local_y = args.index_range_y + (np.asarray(projected_y, dtype=np.float64) - args.min_y)
    projection_time_s = now_s() - projection_start

    assembly_start = now_s()
    targets = np.empty((batch_size, 3), dtype=np.float64)
    targets[:, 0] = local_x
    targets[:, 1] = local_y
    targets[:, 2] = zs
    target_assembly_time_s = now_s() - assembly_start

    timings = {
        "random_generation_time_s": random_generation_time_s,
        "projection_time_s": projection_time_s,
        "target_assembly_time_s": target_assembly_time_s,
    }
    return np.ascontiguousarray(targets, dtype=np.float64), timings


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Random-point BVH-RT benchmark. Core timing excludes random sampling and projection."
    )
    parser.add_argument(
        "--sample-counts",
        nargs="+",
        type=positive_int,
        default=[1_000_000, 10_000_000, 100_000_000],
        help="Point counts to benchmark.",
    )
    parser.add_argument(
        "--query-batch-size",
        type=positive_int,
        default=1_000_000,
        help="Random points generated and queried per outer batch.",
    )
    parser.add_argument(
        "--ray-chunk-size",
        type=nonnegative_int,
        default=0,
        help="Open3D rays per internal cast_rays call; 0 means use --query-batch-size.",
    )
    parser.add_argument("--repeats", type=positive_int, default=1)
    parser.add_argument("--seed", type=int, default=20260517)
    parser.add_argument(
        "--warmup-batches",
        type=nonnegative_int,
        default=0,
        help="Run this many untimed random query batches after scene build and before measured benchmarks.",
    )

    parser.add_argument("--min-lon", type=float, default=DEFAULT_MIN_LON)
    parser.add_argument("--max-lon", type=float, default=DEFAULT_MAX_LON)
    parser.add_argument("--min-lat", type=float, default=DEFAULT_MIN_LAT)
    parser.add_argument("--max-lat", type=float, default=DEFAULT_MAX_LAT)
    parser.add_argument("--fixed-height", type=float, default=30.0)
    parser.add_argument("--random-alt", action="store_true", help="Sample z uniformly from [--min-alt, --max-alt].")
    parser.add_argument("--min-alt", type=float, default=0.0)
    parser.add_argument("--max-alt", type=float, default=100.0)

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
    parser.add_argument("--build-threads", type=int, default=0)
    parser.add_argument("--query-threads", type=int, default=8)
    parser.add_argument("--mesh-offset", nargs=3, type=float, help="Manual offset subtracted from mesh/points.")
    parser.add_argument("--zero-mesh-offset", action="store_true", help="Do not apply Open3D numeric-stability offset.")
    parser.add_argument("--summary-output", default="out/bvh_rt_random_benchmark_summary.csv")
    parser.add_argument("--mesh-report", default="out/bvh_rt_random_benchmark_mesh_report.csv")
    parser.add_argument("--progress-every", type=int, default=0)
    parser.add_argument(
        "--no-batch-progress",
        action="store_true",
        help="Do not print completed outer batches.",
    )
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
    if args.min_lon >= args.max_lon:
        raise ValueError("--min-lon must be less than --max-lon")
    if args.min_lat >= args.max_lat:
        raise ValueError("--min-lat must be less than --max-lat")
    if args.min_alt > args.max_alt:
        raise ValueError("--min-alt must be less than or equal to --max-alt")
    if not args.random_alt and not math.isfinite(args.fixed_height):
        raise ValueError("--fixed-height must be finite")
    if args.query_threads < 0 or args.build_threads < 0:
        raise ValueError("--query-threads and --build-threads must be non-negative")


def finite_or_nan(value: object, default: float = float("nan")) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def run_random_benchmark(
    scene,
    radar_local: np.ndarray,
    mesh_offset: np.ndarray,
    transformer,
    args: argparse.Namespace,
    sample_count: int,
    repeat_index: int,
) -> dict[str, object]:
    seed = args.seed + repeat_index
    rng = np.random.default_rng(seed)
    ray_chunk_size = args.ray_chunk_size if args.ray_chunk_size > 0 else args.query_batch_size

    benchmark_start = now_s()
    processed = 0
    batch_index = 0

    random_generation_time_s = 0.0
    projection_time_s = 0.0
    target_assembly_time_s = 0.0
    core_compute_time_s = 0.0
    cast_rays_time_s = 0.0
    ray_prepare_and_classify_time_s = 0.0
    cast_rays_time_known = True
    prep_time_known = True

    num_visible = 0
    num_blocked = 0
    num_invalid_targets = 0
    num_near_tx_targets = 0
    peak_rss_gb = current_rss_gb()

    while processed < sample_count:
        current_batch_size = min(args.query_batch_size, sample_count - processed)
        targets, sample_timings = generate_random_target_batch(
            rng=rng,
            batch_size=current_batch_size,
            transformer=transformer,
            args=args,
        )
        random_generation_time_s += sample_timings["random_generation_time_s"]
        projection_time_s += sample_timings["projection_time_s"]
        target_assembly_time_s += sample_timings["target_assembly_time_s"]

        compute_start = now_s()
        visible, query_stats = cast_visibility_batch(
            scene=scene,
            tx=radar_local,
            targets=targets,
            offset=mesh_offset,
            eps=args.eps,
            chunk_size=ray_chunk_size,
            query_threads=args.query_threads,
            progress_every=args.progress_every,
        )
        measured_compute_time_s = now_s() - compute_start

        batch_query_time_s = finite_or_nan(query_stats.get("query_time_s"), measured_compute_time_s)
        core_compute_time_s += batch_query_time_s

        batch_cast_time_s = finite_or_nan(query_stats.get("cast_rays_time_s"))
        if math.isfinite(batch_cast_time_s):
            cast_rays_time_s += batch_cast_time_s
        else:
            cast_rays_time_known = False

        batch_prep_time_s = finite_or_nan(query_stats.get("ray_prepare_and_classify_time_s"))
        if math.isfinite(batch_prep_time_s):
            ray_prepare_and_classify_time_s += batch_prep_time_s
        else:
            prep_time_known = False

        num_visible += int(query_stats.get("num_visible", int(visible.sum())))
        num_blocked += int(query_stats.get("num_blocked", current_batch_size - int(visible.sum())))
        num_invalid_targets += int(query_stats.get("num_invalid_targets", 0))
        num_near_tx_targets += int(query_stats.get("num_near_tx_targets", 0))
        peak_rss_gb = max(peak_rss_gb, finite_or_nan(query_stats.get("peak_rss_gb"), current_rss_gb()))

        processed += current_batch_size
        batch_index += 1
        if not args.no_batch_progress:
            print(
                f"  sample_count={sample_count} repeat={repeat_index + 1}/{args.repeats} "
                f"batch={batch_index} processed={processed}/{sample_count} "
                f"core_compute_time_s={core_compute_time_s:.6f}"
            )

        del targets, visible

    benchmark_wall_time_s = now_s() - benchmark_start
    throughput_qps = float(sample_count / core_compute_time_s) if core_compute_time_s > 0.0 else float("nan")
    cast_rays_throughput_qps = (
        float(sample_count / cast_rays_time_s)
        if cast_rays_time_known and cast_rays_time_s > 0.0
        else float("nan")
    )
    visible_ratio = float(num_visible / sample_count) if sample_count else float("nan")
    sampling_excluded_time_s = random_generation_time_s + projection_time_s + target_assembly_time_s

    return {
        "method": "BVH-RT-Open3D-random",
        "sample_count": sample_count,
        "repeat_index": repeat_index,
        "seed": seed,
        "min_lon": args.min_lon,
        "max_lon": args.max_lon,
        "min_lat": args.min_lat,
        "max_lat": args.max_lat,
        "random_alt": int(args.random_alt),
        "fixed_height": args.fixed_height if not args.random_alt else float("nan"),
        "min_alt": args.min_alt,
        "max_alt": args.max_alt,
        "radar_lon": args.radar_lon,
        "radar_lat": args.radar_lat,
        "radar_alt": args.radar_alt,
        "radar_local_x": radar_local[0],
        "radar_local_y": radar_local[1],
        "radar_local_z": radar_local[2],
        "query_batch_size": args.query_batch_size,
        "ray_chunk_size": ray_chunk_size,
        "build_threads": args.build_threads,
        "query_threads": args.query_threads,
        "eps": args.eps,
        "random_generation_time_s": random_generation_time_s,
        "projection_time_s": projection_time_s,
        "target_assembly_time_s": target_assembly_time_s,
        "sampling_excluded_time_s": sampling_excluded_time_s,
        "core_compute_time_s": core_compute_time_s,
        "query_time_s": core_compute_time_s,
        "cast_rays_time_s": cast_rays_time_s if cast_rays_time_known else float("nan"),
        "ray_prepare_and_classify_time_s": (
            ray_prepare_and_classify_time_s if prep_time_known else float("nan")
        ),
        "throughput_qps": throughput_qps,
        "cast_rays_throughput_qps": cast_rays_throughput_qps,
        "num_visible": num_visible,
        "num_blocked": num_blocked,
        "num_invalid_targets": num_invalid_targets,
        "num_near_tx_targets": num_near_tx_targets,
        "visible_ratio": visible_ratio,
        "peak_rss_gb": peak_rss_gb,
        "benchmark_wall_time_s": benchmark_wall_time_s,
    }


def run_warmup_batches(
    scene,
    radar_local: np.ndarray,
    mesh_offset: np.ndarray,
    transformer,
    args: argparse.Namespace,
) -> None:
    if args.warmup_batches <= 0:
        return

    rng = np.random.default_rng(args.seed - 1)
    ray_chunk_size = args.ray_chunk_size if args.ray_chunk_size > 0 else args.query_batch_size
    print(f"Running {args.warmup_batches} untimed warm-up batch(es)...")
    for batch_index in range(args.warmup_batches):
        targets, _ = generate_random_target_batch(
            rng=rng,
            batch_size=args.query_batch_size,
            transformer=transformer,
            args=args,
        )
        warmup_start = now_s()
        visible, _ = cast_visibility_batch(
            scene=scene,
            tx=radar_local,
            targets=targets,
            offset=mesh_offset,
            eps=args.eps,
            chunk_size=ray_chunk_size,
            query_threads=args.query_threads,
            progress_every=0,
        )
        print(
            f"  warmup batch {batch_index + 1}/{args.warmup_batches}: "
            f"{now_s() - warmup_start:.6f} s"
        )
        del targets, visible


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    validate_args(args)

    Transformer = import_pyproj()
    transformer = Transformer.from_crs(args.source_crs, args.projected_crs, always_xy=True)

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

    common_scene_row = {
        "num_obj_files": build_stats["num_obj_files"],
        "total_vertices": build_stats["total_vertices"],
        "total_triangles": build_stats["total_triangles"],
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
        "bbox_scan_time_s": bbox_scan_time_s,
        "mesh_read_time_s": build_stats["mesh_read_time_s"],
        "mesh_clean_shift_time_s": build_stats["mesh_clean_shift_time_s"],
        "scene_add_time_s": build_stats["scene_add_time_s"],
        "scene_build_time_s": build_stats["scene_build_time_s"],
        "rss_before_build_gb": rss_before_build_gb,
        "rss_after_build_gb": rss_after_build_gb,
    }

    run_warmup_batches(
        scene=scene,
        radar_local=radar_local,
        mesh_offset=mesh_offset,
        transformer=transformer,
        args=args,
    )

    rows: list[dict[str, object]] = []
    summary_path = Path(args.summary_output)
    for sample_count in args.sample_counts:
        for repeat_index in range(args.repeats):
            print(
                f"Running random BVH benchmark: sample_count={sample_count}, "
                f"repeat={repeat_index + 1}/{args.repeats}, query_threads={args.query_threads}"
            )
            row = run_random_benchmark(
                scene=scene,
                radar_local=radar_local,
                mesh_offset=mesh_offset,
                transformer=transformer,
                args=args,
                sample_count=sample_count,
                repeat_index=repeat_index,
            )
            row.update(common_scene_row)
            rows.append(row)
            write_csv(summary_path, rows)

            print(
                f"  done sample_count={sample_count}: "
                f"core_compute_time_s={row['core_compute_time_s']:.6f}, "
                f"cast_rays_time_s={row['cast_rays_time_s']:.6f}, "
                f"throughput_qps={row['throughput_qps']:.2f}, "
                f"visible_ratio={row['visible_ratio']:.6f}"
            )

    meta_path = summary_path.with_suffix(summary_path.suffix + ".json")
    with meta_path.open("w", encoding="utf-8") as f:
        json.dump(rows, f, ensure_ascii=False, indent=2)

    print(f"Summary written to {summary_path}")
    print(f"Summary JSON written to {meta_path}")
    print(f"Mesh report written to {args.mesh_report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
