#!/usr/bin/env python3
"""BVH-RT multi-OBJ 3D visibility experiments.

The default command runs an Open3D RaycastingScene benchmark. Extra
subcommands generate reproducible target-point files and compare visibility
label arrays. All coordinates used by ``run`` are expected to be in the same
local Cartesian frame as the OBJ vertices.
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
from typing import Iterable, Sequence

import numpy as np

try:
    from PIL import Image
except ImportError:  # pragma: no cover - optional dependency
    Image = None

try:
    import psutil
except ImportError:  # pragma: no cover - optional dependency
    psutil = None


METHOD_NAME = "BVH-RT-Open3D"
DEFAULT_SIZES = "100000,500000,1000000,5000000,10000000"


def now_s() -> float:
    return time.perf_counter()


def current_rss_gb() -> float:
    if psutil is None:
        return float("nan")
    return psutil.Process().memory_info().rss / (1024.0**3)


def import_open3d():
    try:
        import open3d as o3d
    except ImportError as exc:  # pragma: no cover - user environment dependent
        raise SystemExit(
            "Open3D is required for the BVH-RT run command. Install it with: "
            "python3.10 -m pip install open3d numpy psutil. "
            "Use Python 3.10/3.11/3.12; Open3D may not provide wheels for newer Python versions."
        ) from exc
    return o3d


_O3D = None


def get_open3d():
    global _O3D
    if _O3D is None:
        _O3D = import_open3d()
    return _O3D


def parse_size_list(text: str) -> list[int]:
    sizes: list[int] = []
    for item in text.split(","):
        item = item.strip().replace("_", "")
        if not item:
            continue
        value = int(item)
        if value <= 0:
            raise ValueError("Target sizes must be positive.")
        sizes.append(value)
    if not sizes:
        raise ValueError("Target size list is empty.")
    return sizes


def parse_float_list(text: str) -> list[float]:
    values = [float(item.strip()) for item in text.split(",") if item.strip()]
    if not values:
        raise ValueError("Float list is empty.")
    return values


def discover_obj_paths(
    obj_args: Sequence[str],
    obj_dir: str | None,
    obj_glob: str,
) -> list[Path]:
    paths: list[Path] = []

    for item in obj_args:
        p = Path(item)
        if any(ch in item for ch in ["*", "?", "["]):
            paths.extend(Path(match) for match in sorted(glob.glob(item, recursive=True)))
        elif p.is_dir():
            paths.extend(sorted(p.rglob(obj_glob)))
        else:
            paths.append(p)

    if obj_dir:
        paths.extend(sorted(Path(obj_dir).rglob(obj_glob)))

    unique: list[Path] = []
    seen: set[str] = set()
    for p in paths:
        if p.suffix.lower() != ".obj":
            continue
        resolved = p.resolve(strict=False)
        key = str(resolved)
        if key not in seen:
            unique.append(resolved)
            seen.add(key)

    missing = [str(p) for p in unique if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing OBJ files:\n" + "\n".join(missing[:20]))
    if not unique:
        raise ValueError("No OBJ files found. Use --obj-dir or --obj.")
    return unique


def scan_obj_bbox(obj_paths: Iterable[Path]) -> tuple[np.ndarray, np.ndarray, int, float]:
    """Scan OBJ vertex records without constructing full meshes."""
    bbox_min = np.array([np.inf, np.inf, np.inf], dtype=np.float64)
    bbox_max = np.array([-np.inf, -np.inf, -np.inf], dtype=np.float64)
    vertex_count = 0
    start = now_s()

    for path in obj_paths:
        with path.open("r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if not line.startswith("v "):
                    continue
                parts = line.split()
                if len(parts) < 4:
                    continue
                try:
                    vertex = np.array(
                        [float(parts[1]), float(parts[2]), float(parts[3])],
                        dtype=np.float64,
                    )
                except ValueError:
                    continue
                bbox_min = np.minimum(bbox_min, vertex)
                bbox_max = np.maximum(bbox_max, vertex)
                vertex_count += 1

    elapsed = now_s() - start
    if vertex_count == 0:
        raise ValueError("No OBJ vertices found while scanning bbox.")
    return bbox_min, bbox_max, vertex_count, elapsed


def read_mesh(o3d, path: Path):
    mesh = o3d.io.read_triangle_mesh(str(path), enable_post_processing=True)
    if len(mesh.vertices) == 0:
        raise ValueError(f"OBJ has no vertices: {path}")
    if len(mesh.triangles) == 0:
        raise ValueError(f"OBJ has no triangles after loading: {path}")
    return mesh


def clean_mesh(mesh):
    mesh.remove_degenerate_triangles()
    mesh.remove_duplicated_triangles()
    mesh.remove_duplicated_vertices()
    mesh.remove_unreferenced_vertices()
    return mesh


def build_scene(
    o3d,
    obj_paths: list[Path],
    offset: np.ndarray,
    build_threads: int,
) -> tuple[object, dict[str, float | int], list[dict[str, object]]]:
    scene = o3d.t.geometry.RaycastingScene(nthreads=build_threads)

    total_vertices = 0
    total_triangles = 0
    read_time_s = 0.0
    clean_time_s = 0.0
    add_time_s = 0.0
    mesh_rows: list[dict[str, object]] = []

    scene_start = now_s()
    for path in obj_paths:
        t0 = now_s()
        mesh = read_mesh(o3d, path)
        read_elapsed = now_s() - t0
        read_time_s += read_elapsed

        t1 = now_s()
        mesh = clean_mesh(mesh)
        vertices = np.asarray(mesh.vertices, dtype=np.float64)
        if not np.isfinite(vertices).all():
            raise ValueError(f"OBJ contains non-finite vertices: {path}")
        mesh.vertices = o3d.utility.Vector3dVector(vertices - offset)
        n_vertices = len(mesh.vertices)
        n_triangles = len(mesh.triangles)
        clean_elapsed = now_s() - t1
        clean_time_s += clean_elapsed

        t2 = now_s()
        tensor_mesh = o3d.t.geometry.TriangleMesh.from_legacy(mesh)
        geom_id = scene.add_triangles(tensor_mesh)
        add_elapsed = now_s() - t2
        add_time_s += add_elapsed

        total_vertices += n_vertices
        total_triangles += n_triangles
        mesh_rows.append(
            {
                "obj_file": str(path),
                "geometry_id": int(geom_id),
                "vertices": n_vertices,
                "triangles": n_triangles,
                "read_time_s": read_elapsed,
                "clean_shift_time_s": clean_elapsed,
                "add_to_scene_time_s": add_elapsed,
            }
        )

    scene_build_time_s = now_s() - scene_start
    stats: dict[str, float | int] = {
        "num_obj_files": len(obj_paths),
        "total_vertices": total_vertices,
        "total_triangles": total_triangles,
        "mesh_read_time_s": read_time_s,
        "mesh_clean_shift_time_s": clean_time_s,
        "scene_add_time_s": add_time_s,
        "scene_build_time_s": scene_build_time_s,
    }
    return scene, stats, mesh_rows


def load_targets(path: Path, mmap_npy: bool) -> np.ndarray:
    suffix = path.suffix.lower()
    if suffix == ".npy":
        arr = np.load(path, mmap_mode="r" if mmap_npy else None, allow_pickle=False)
    elif suffix == ".npz":
        data = np.load(path, allow_pickle=False)
        key = "targets" if "targets" in data else sorted(data.files)[0]
        arr = data[key]
    elif suffix in {".csv", ".txt"}:
        try:
            arr = np.loadtxt(path, delimiter=",", comments="#", dtype=np.float64)
        except ValueError:
            named = np.genfromtxt(path, delimiter=",", names=True, dtype=np.float64)
            cols = named.dtype.names
            if cols is None or len(cols) < 3:
                raise ValueError(f"Target CSV must contain at least x,y,z columns: {path}")
            arr = np.vstack([named[cols[0]], named[cols[1]], named[cols[2]]]).T
    else:
        raise ValueError(f"Unsupported target file format: {path}")

    arr = np.asarray(arr)
    if arr.ndim == 1 and arr.shape[0] >= 3:
        arr = arr.reshape(1, arr.shape[0])
    if arr.ndim != 2 or arr.shape[1] < 3:
        raise ValueError(f"Target array must have shape (N,3) or more columns: {path}")
    return arr


def load_points(path: Path, mmap_npy: bool = True) -> np.ndarray:
    return load_targets(path, mmap_npy=mmap_npy)


def save_visibility_result(
    output_path: Path,
    targets: np.ndarray,
    visible: np.ndarray,
    method: str,
    target_file: Path,
    summary: dict[str, object],
    include_points: bool,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    visible_arr = np.asarray(visible, dtype=np.bool_)
    if include_points:
        points = np.asarray(targets[:, 0:3], dtype=np.float64)
        np.savez_compressed(
            output_path,
            points=points,
            visible=visible_arr,
            blocked=~visible_arr,
            method=np.array(method),
            target_file=np.array(str(target_file)),
            summary_json=np.array(json.dumps(summary, ensure_ascii=False)),
        )
    else:
        np.savez_compressed(
            output_path,
            visible=visible_arr,
            blocked=~visible_arr,
            method=np.array(method),
            target_file=np.array(str(target_file)),
            summary_json=np.array(json.dumps(summary, ensure_ascii=False)),
        )


def read_summary_json_from_npz(path: Path) -> dict[str, object]:
    if path.suffix.lower() != ".npz":
        return {}
    data = np.load(path, allow_pickle=False)
    if "summary_json" not in data.files:
        return {}
    raw = str(np.asarray(data["summary_json"]).item())
    try:
        parsed = json.loads(raw)
    except json.JSONDecodeError:
        return {}
    return parsed if isinstance(parsed, dict) else {}


def cast_visibility_batch(
    scene,
    tx: np.ndarray,
    targets: np.ndarray,
    offset: np.ndarray,
    eps: float,
    chunk_size: int,
    query_threads: int,
    progress_every: int,
) -> tuple[np.ndarray, dict[str, float | int]]:
    if chunk_size <= 0:
        raise ValueError("--chunk-size must be positive.")
    if eps < 0:
        raise ValueError("--eps must be non-negative.")

    tx_local = np.asarray(tx, dtype=np.float64) - offset
    if tx_local.shape != (3,) or not np.isfinite(tx_local).all():
        raise ValueError("Transmitter point must contain three finite coordinates.")

    n = int(targets.shape[0])
    visible = np.empty(n, dtype=np.bool_)
    peak_rss_gb = current_rss_gb()
    invalid_targets = 0
    near_tx_targets = 0
    cast_rays_time_s = 0.0

    start_time = now_s()
    for start in range(0, n, chunk_size):
        end = min(start + chunk_size, n)
        pts = np.asarray(targets[start:end, 0:3], dtype=np.float64)
        pts_local = pts - offset

        finite = np.isfinite(pts_local).all(axis=1)
        vec = pts_local - tx_local[None, :]
        dist = np.linalg.norm(vec, axis=1)
        valid = finite & (dist > (2.0 * eps))
        near_tx = finite & ~valid

        dirs = np.zeros_like(vec)
        dirs[valid] = vec[valid] / dist[valid, None]
        dirs[~valid] = np.array([1.0, 0.0, 0.0], dtype=np.float64)
        origins = tx_local[None, :] + eps * dirs

        rays = np.empty((end - start, 6), dtype=np.float32)
        rays[:, 0:3] = origins.astype(np.float32)
        rays[:, 3:6] = dirs.astype(np.float32)

        cast_start = now_s()
        ans = scene.cast_rays(get_open3d().core.Tensor(rays), nthreads=query_threads)
        cast_rays_time_s += now_s() - cast_start
        t_hit = ans["t_hit"].numpy()

        blocked = np.isfinite(t_hit) & (t_hit < (dist - 2.0 * eps))
        blocked[~valid] = False
        chunk_visible = ~blocked
        chunk_visible[~finite] = False
        chunk_visible[near_tx] = True
        visible[start:end] = chunk_visible

        invalid_targets += int((~finite).sum())
        near_tx_targets += int(near_tx.sum())

        rss = current_rss_gb()
        if not math.isnan(rss):
            peak_rss_gb = max(peak_rss_gb, rss)
        if progress_every > 0 and (end == n or end // progress_every > start // progress_every):
            print(f"  cast {end}/{n} rays")

    query_time_s = now_s() - start_time
    num_visible = int(visible.sum())
    num_blocked = int(n - num_visible - invalid_targets)

    stats: dict[str, float | int] = {
        "num_queries": n,
        "query_time_s": query_time_s,
        "cast_rays_time_s": cast_rays_time_s,
        "ray_prepare_and_classify_time_s": max(0.0, query_time_s - cast_rays_time_s),
        "throughput_qps": float(n / query_time_s) if query_time_s > 0 else float("inf"),
        "cast_rays_throughput_qps": float(n / cast_rays_time_s) if cast_rays_time_s > 0 else float("inf"),
        "num_visible": num_visible,
        "num_blocked": num_blocked,
        "num_invalid_targets": invalid_targets,
        "num_near_tx_targets": near_tx_targets,
        "visible_ratio": float(num_visible / n) if n else float("nan"),
        "peak_rss_gb": peak_rss_gb,
    }
    return visible, stats


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def append_csv_row(path: Path, row: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    exists = path.exists()
    with path.open("a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(row.keys()))
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def build_height_values(args: argparse.Namespace) -> np.ndarray:
    if args.height_list:
        heights = np.asarray(parse_float_list(args.height_list), dtype=np.float64)
    else:
        start, stop, step = args.height_range
        if step <= 0:
            raise ValueError("--height-range STEP must be positive.")
        heights = np.arange(start, stop, step, dtype=np.float64)
    if heights.size == 0:
        raise ValueError("No heights were generated.")
    if not np.isfinite(heights).all():
        raise ValueError("Heights must be finite.")
    return heights


def fill_points_from_xy_indices(
    out: np.ndarray,
    xs: np.ndarray,
    ys: np.ndarray,
    xy_indices: np.ndarray,
    z_value: float,
    offset: int = 0,
) -> None:
    nx = xs.size
    out_slice = out[offset : offset + xy_indices.size]
    out_slice[:, 0] = xs[xy_indices % nx]
    out_slice[:, 1] = ys[xy_indices // nx]
    out_slice[:, 2] = z_value


def fill_points_from_global_range(
    out: np.ndarray,
    xs: np.ndarray,
    ys: np.ndarray,
    heights: np.ndarray,
    total_size: int,
    chunk_rows: int,
) -> None:
    nx = xs.size
    nxy = xs.size * ys.size
    for start in range(0, total_size, chunk_rows):
        end = min(start + chunk_rows, total_size)
        global_idx = np.arange(start, end, dtype=np.int64)
        z_idx = global_idx // nxy
        xy_idx = global_idx % nxy
        out_chunk = out[start:end]
        out_chunk[:, 0] = xs[xy_idx % nx]
        out_chunk[:, 1] = ys[xy_idx // nx]
        out_chunk[:, 2] = heights[z_idx]


def generate_targets(args: argparse.Namespace) -> int:
    if args.bbox:
        x_min, x_max, y_min, y_max = args.bbox
        bbox_source = "manual"
    else:
        obj_paths = discover_obj_paths(args.obj, args.obj_dir, args.obj_glob)
        bbox_min, bbox_max, _, _ = scan_obj_bbox(obj_paths)
        x_min, y_min = float(bbox_min[0]), float(bbox_min[1])
        x_max, y_max = float(bbox_max[0]), float(bbox_max[1])
        bbox_source = "obj"

    if args.bbox_padding != 0:
        x_min -= args.bbox_padding
        y_min -= args.bbox_padding
        x_max += args.bbox_padding
        y_max += args.bbox_padding
    if x_min >= x_max or y_min >= y_max:
        raise ValueError("Invalid bbox: min values must be less than max values.")
    if args.spacing <= 0:
        raise ValueError("--spacing must be positive.")

    xs = np.arange(x_min, x_max, args.spacing, dtype=np.float64)
    ys = np.arange(y_min, y_max, args.spacing, dtype=np.float64)
    if xs.size == 0 or ys.size == 0:
        raise ValueError("The bbox and spacing generate an empty grid.")

    heights = build_height_values(args)
    nxy = int(xs.size * ys.size)
    total_available = int(nxy * heights.size)
    sizes = parse_size_list(args.sizes)
    if max(sizes) > total_available:
        raise ValueError(
            f"Requested {max(sizes)} points but only {total_available} are available "
            f"({xs.size} x {ys.size} x {heights.size})."
        )

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    manifest_rows: list[dict[str, object]] = []
    dtype = np.dtype(args.dtype)

    for size in sizes:
        path = out_dir / f"{args.prefix}_{size}.npy"
        if path.exists() and not args.overwrite:
            raise FileExistsError(f"Target file exists, pass --overwrite to replace: {path}")

        out = np.lib.format.open_memmap(path, mode="w+", dtype=dtype, shape=(size, 3))
        sampled = False
        if args.sample_smallest and size == min(sizes) and size <= nxy:
            xy_idx = rng.choice(nxy, size=size, replace=False)
            fill_points_from_xy_indices(out, xs, ys, xy_idx.astype(np.int64), float(heights[0]))
            sampled = True
        else:
            fill_points_from_global_range(out, xs, ys, heights, size, args.chunk_rows)
        out.flush()

        row = {
            "target_file": str(path),
            "num_points": size,
            "dtype": str(dtype),
            "bbox_source": bbox_source,
            "x_min": x_min,
            "x_max": x_max,
            "y_min": y_min,
            "y_max": y_max,
            "spacing": args.spacing,
            "height_min": float(heights.min()),
            "height_max": float(heights.max()),
            "num_heights": int(heights.size),
            "sampled": sampled,
            "seed": args.seed,
        }
        manifest_rows.append(row)
        print(f"Wrote {size} targets: {path}")

    if args.manifest:
        write_csv(Path(args.manifest), manifest_rows)
        print(f"Manifest written to {args.manifest}")
    return 0


def generate_candidates(args: argparse.Namespace) -> int:
    x_min, x_max, y_min, y_max = args.bbox
    if x_min >= x_max or y_min >= y_max:
        raise ValueError("Invalid bbox: min values must be less than max values.")
    if args.spacing <= 0:
        raise ValueError("--spacing must be positive.")

    xs = np.arange(x_min, x_max, args.spacing, dtype=np.float64)
    ys = np.arange(y_min, y_max, args.spacing, dtype=np.float64)
    if xs.size == 0 or ys.size == 0:
        raise ValueError("The bbox and spacing generate an empty candidate grid.")

    total = int(xs.size * ys.size)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    dtype = np.dtype(args.dtype)
    points = np.lib.format.open_memmap(output, mode="w+", dtype=dtype, shape=(total, 3))
    idx = 0
    for y in ys:
        end = idx + xs.size
        points[idx:end, 0] = xs
        points[idx:end, 1] = y
        points[idx:end, 2] = args.z
        idx = end
    points.flush()

    row = {
        "candidate_file": str(output),
        "num_points": total,
        "width": int(xs.size),
        "height": int(ys.size),
        "x_min": x_min,
        "x_max": x_max,
        "y_min": y_min,
        "y_max": y_max,
        "z": args.z,
        "spacing": args.spacing,
        "dtype": str(dtype),
    }
    if args.manifest:
        append_csv_row(Path(args.manifest), row)
    print(json.dumps(row, ensure_ascii=False, indent=2))
    return 0


def run_bvh_rt(args: argparse.Namespace) -> int:
    o3d = get_open3d()
    obj_paths = discover_obj_paths(args.obj, args.obj_dir, args.obj_glob)
    print(f"Discovered {len(obj_paths)} OBJ files.")

    bbox_scan_time_s = 0.0
    bbox_vertex_count = 0
    bbox_min = np.array([float("nan")] * 3, dtype=np.float64)
    bbox_max = np.array([float("nan")] * 3, dtype=np.float64)
    if args.offset is not None:
        offset = np.array(args.offset, dtype=np.float64)
    elif args.zero_offset:
        offset = np.zeros(3, dtype=np.float64)
    else:
        bbox_min, bbox_max, bbox_vertex_count, bbox_scan_time_s = scan_obj_bbox(obj_paths)
        offset = (bbox_min + bbox_max) * 0.5

    if offset.shape != (3,) or not np.isfinite(offset).all():
        raise ValueError("Coordinate offset must contain three finite values.")
    print(f"Using coordinate offset: {offset.tolist()}")

    rss_before_build_gb = current_rss_gb()
    scene, build_stats, mesh_rows = build_scene(o3d, obj_paths, offset, args.build_threads)
    rss_after_build_gb = current_rss_gb()
    write_csv(Path(args.mesh_report), mesh_rows)
    print(f"Scene built with {build_stats['total_triangles']} triangles.")
    print(f"Mesh report written to {args.mesh_report}")

    tx = np.array(args.tx, dtype=np.float64)
    output_path = Path(args.output)
    labels_dir = Path(args.labels_dir) if args.labels_dir else None
    if labels_dir:
        labels_dir.mkdir(parents=True, exist_ok=True)

    target_paths = [Path(p).resolve(strict=False) for p in args.targets]
    for target_path in target_paths:
        load_start = now_s()
        targets = load_targets(target_path, mmap_npy=not args.no_mmap)
        target_load_time_s = now_s() - load_start
        print(f"Loaded {targets.shape[0]} targets from {target_path}")

        for rep in range(args.repeat):
            visible, query_stats = cast_visibility_batch(
                scene=scene,
                tx=tx,
                targets=targets,
                offset=offset,
                eps=args.eps,
                chunk_size=args.chunk_size,
                query_threads=args.query_threads,
                progress_every=args.progress_every,
            )

            if labels_dir:
                label_path = labels_dir / f"{target_path.stem}_rep{rep + 1}_visible.npy"
                np.save(label_path, visible)
            else:
                label_path = ""

            row: dict[str, object] = {
                "method": METHOD_NAME,
                "target_file": str(target_path),
                "repeat": rep + 1,
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
                "offset_x": offset[0],
                "offset_y": offset[1],
                "offset_z": offset[2],
                "tx_x": tx[0],
                "tx_y": tx[1],
                "tx_z": tx[2],
                "eps": args.eps,
                "chunk_size": args.chunk_size,
                "build_threads": args.build_threads,
                "query_threads": args.query_threads,
                "bbox_scan_time_s": bbox_scan_time_s,
                "mesh_read_time_s": build_stats["mesh_read_time_s"],
                "mesh_clean_shift_time_s": build_stats["mesh_clean_shift_time_s"],
                "scene_add_time_s": build_stats["scene_add_time_s"],
                "scene_build_time_s": build_stats["scene_build_time_s"],
                "rss_before_build_gb": rss_before_build_gb,
                "rss_after_build_gb": rss_after_build_gb,
                "target_load_time_s": target_load_time_s,
                "num_queries": query_stats["num_queries"],
                "query_time_s": query_stats["query_time_s"],
                "throughput_qps": query_stats["throughput_qps"],
                "num_visible": query_stats["num_visible"],
                "num_blocked": query_stats["num_blocked"],
                "num_invalid_targets": query_stats["num_invalid_targets"],
                "num_near_tx_targets": query_stats["num_near_tx_targets"],
                "visible_ratio": query_stats["visible_ratio"],
                "peak_rss_gb": query_stats["peak_rss_gb"],
                "labels_file": str(label_path),
                "result_file": "",
            }
            if args.results_dir:
                results_dir = Path(args.results_dir)
                results_dir.mkdir(parents=True, exist_ok=True)
                result_path = results_dir / f"{target_path.stem}_rep{rep + 1}_{args.result_suffix}.npz"
                save_visibility_result(
                    output_path=result_path,
                    targets=targets,
                    visible=visible,
                    method=METHOD_NAME,
                    target_file=target_path,
                    summary=row,
                    include_points=args.result_include_points,
                )
                row["result_file"] = str(result_path)
            append_csv_row(output_path, row)
            print(
                f"{target_path.name} rep {rep + 1}: "
                f"{query_stats['query_time_s']:.6f}s, "
                f"{query_stats['throughput_qps']:.2f} q/s, "
                f"visible={query_stats['num_visible']}, "
                f"blocked={query_stats['num_blocked']}"
            )

    print(f"Summary written to {output_path}")
    return 0


def load_label_array(path: Path) -> np.ndarray:
    suffix = path.suffix.lower()
    if suffix == ".npy":
        arr = np.load(path, allow_pickle=False)
    elif suffix == ".npz":
        data = np.load(path, allow_pickle=False)
        key = "visible" if "visible" in data else sorted(data.files)[0]
        arr = data[key]
    elif suffix in {".csv", ".txt"}:
        arr = np.loadtxt(path, delimiter=",", comments="#")
    else:
        raise ValueError(f"Unsupported label file format: {path}")
    arr = np.asarray(arr)
    if arr.ndim != 1:
        arr = arr.reshape(-1)
    if arr.dtype == np.bool_:
        return arr
    return arr > 0


def load_result_or_labels(path: Path) -> tuple[np.ndarray, dict[str, object]]:
    suffix = path.suffix.lower()
    if suffix == ".npz":
        data = np.load(path, allow_pickle=False)
        if "visible" in data.files:
            labels = np.asarray(data["visible"], dtype=np.bool_).reshape(-1)
            return labels, read_summary_json_from_npz(path)
    return load_label_array(path), {}


def local_to_lonlat_linear(points: np.ndarray, args: argparse.Namespace) -> np.ndarray:
    pts = np.asarray(points[:, 0:3], dtype=np.float64)
    lon = args.origin_lon + (pts[:, 0] - args.origin_x) * args.meters_to_lon
    lat = args.origin_lat + (pts[:, 1] - args.origin_y) * args.meters_to_lat
    return np.column_stack([lon, lat, pts[:, 2]])


def local_to_lonlat_pyproj(points: np.ndarray, args: argparse.Namespace) -> np.ndarray:
    try:
        from pyproj import Transformer
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise SystemExit("pyproj is required for --mode pyproj. Install it with: pip install pyproj") from exc

    pts = np.asarray(points[:, 0:3], dtype=np.float64)
    x = pts[:, 0] + args.projected_offset_x
    y = pts[:, 1] + args.projected_offset_y
    transformer = Transformer.from_crs(args.source_crs, args.target_crs, always_xy=True)
    lon, lat = transformer.transform(x, y)
    return np.column_stack([lon, lat, pts[:, 2]])


def convert_local_to_lonlat(args: argparse.Namespace) -> int:
    points = load_points(Path(args.points), mmap_npy=not args.no_mmap)
    if args.mode == "linear":
        lonlat = local_to_lonlat_linear(points, args)
    else:
        lonlat = local_to_lonlat_pyproj(points, args)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.suffix.lower() == ".npy":
        np.save(output, lonlat.astype(np.float64))
    elif output.suffix.lower() == ".npz":
        np.savez_compressed(output, lonlat=lonlat.astype(np.float64), source_points=str(Path(args.points)))
    else:
        header = "lon,lat,z"
        np.savetxt(output, lonlat, delimiter=",", header=header, comments="", fmt="%.10f")
    print(f"Wrote converted lon/lat points to {output}")
    return 0


def merge_visibility_result(args: argparse.Namespace) -> int:
    points = load_points(Path(args.points), mmap_npy=not args.no_mmap)
    labels = load_label_array(Path(args.labels))
    if points.shape[0] != labels.shape[0]:
        raise ValueError(f"Point/label count mismatch: {points.shape[0]} vs {labels.shape[0]}")

    visible = labels.astype(np.bool_)
    n = int(visible.size)
    summary = {
        "method": args.method,
        "target_file": str(Path(args.points).resolve(strict=False)),
        "label_file": str(Path(args.labels).resolve(strict=False)),
        "num_queries": n,
        "num_visible": int(visible.sum()),
        "num_blocked": int(n - visible.sum()),
        "visible_ratio": float(visible.mean()) if n else float("nan"),
        "query_time_s": args.query_time_s,
        "throughput_qps": float(n / args.query_time_s) if args.query_time_s and args.query_time_s > 0 else float("nan"),
        "peak_rss_gb": args.peak_rss_gb,
        "threads": args.threads,
    }
    save_visibility_result(
        output_path=Path(args.output),
        targets=points,
        visible=visible,
        method=args.method,
        target_file=Path(args.points),
        summary=summary,
        include_points=args.include_points,
    )
    if args.summary_csv:
        row = dict(summary)
        row["result_file"] = str(Path(args.output))
        append_csv_row(Path(args.summary_csv), row)
    print(f"Wrote visibility result to {args.output}")
    return 0


def write_png(path: Path, rgb: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if Image is not None:
        Image.fromarray(rgb, mode="RGB").save(path)
        return
    try:
        import imageio.v2 as imageio
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise SystemExit("Install pillow or imageio to write PNG files: pip install pillow") from exc
    imageio.imwrite(path, rgb)


def normalize_grid_axis(values: np.ndarray, spacing: float | None, name: str) -> tuple[np.ndarray, np.ndarray]:
    rounded = np.round(values.astype(np.float64), decimals=9)
    unique = np.unique(rounded)
    if unique.size == 0:
        raise ValueError(f"No {name} coordinates found.")
    if spacing is not None and spacing > 0:
        start = unique.min()
        stop = unique.max() + spacing * 0.5
        axis = np.round(np.arange(start, stop, spacing, dtype=np.float64), decimals=9)
    else:
        axis = unique
    indices = np.searchsorted(axis, rounded)
    if np.any(indices < 0) or np.any(indices >= axis.size) or np.any(axis[indices] != rounded):
        raise ValueError(f"{name} coordinates do not fit the inferred regular grid.")
    return axis, indices


def visualize_slice(args: argparse.Namespace) -> int:
    points = load_points(Path(args.points), mmap_npy=not args.no_mmap)
    labels, _ = load_result_or_labels(Path(args.result))
    if points.shape[0] != labels.shape[0]:
        raise ValueError(f"Point/result count mismatch: {points.shape[0]} vs {labels.shape[0]}")

    pts = np.asarray(points[:, 0:3], dtype=np.float64)
    if args.z is None:
        z_values = np.unique(np.round(pts[:, 2], decimals=9))
        if z_values.size != 1:
            raise ValueError("Multiple z slices found. Pass --z to select one slice.")
        mask = np.ones(pts.shape[0], dtype=np.bool_)
        selected_z = float(z_values[0])
    else:
        selected_z = args.z
        mask = np.abs(pts[:, 2] - selected_z) <= args.z_tol
    if not np.any(mask):
        raise ValueError(f"No points found at z={selected_z} within tolerance {args.z_tol}")

    slice_pts = pts[mask]
    slice_labels = labels[mask]
    xs, x_idx = normalize_grid_axis(slice_pts[:, 0], args.spacing, "x")
    ys, y_idx = normalize_grid_axis(slice_pts[:, 1], args.spacing, "y")
    rgb = np.full((ys.size, xs.size, 3), 230, dtype=np.uint8)

    visible_color = np.array(args.visible_color, dtype=np.uint8)
    blocked_color = np.array(args.blocked_color, dtype=np.uint8)
    rgb[ys.size - 1 - y_idx, x_idx] = np.where(slice_labels[:, None], visible_color, blocked_color)
    write_png(Path(args.output), rgb)

    summary = {
        "points": int(slice_labels.size),
        "z": selected_z,
        "width": int(xs.size),
        "height": int(ys.size),
        "visible": int(slice_labels.sum()),
        "blocked": int(slice_labels.size - slice_labels.sum()),
        "visible_ratio": float(slice_labels.mean()),
    }
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


def compare_results(args: argparse.Namespace) -> int:
    a_labels, a_summary = load_result_or_labels(Path(args.a))
    b_labels, b_summary = load_result_or_labels(Path(args.b))
    if a_labels.shape != b_labels.shape:
        raise ValueError(f"Result shapes differ: {a_labels.shape} vs {b_labels.shape}")

    n = int(a_labels.size)
    both_visible = int(np.count_nonzero(a_labels & b_labels))
    both_blocked = int(np.count_nonzero(~a_labels & ~b_labels))
    a_only = int(np.count_nonzero(a_labels & ~b_labels))
    b_only = int(np.count_nonzero(~a_labels & b_labels))
    agreement = float((both_visible + both_blocked) / n) if n else float("nan")

    a_time = float(a_summary.get("query_time_s", float("nan")))
    b_time = float(b_summary.get("query_time_s", float("nan")))
    row = {
        "file_a": str(Path(args.a).resolve(strict=False)),
        "file_b": str(Path(args.b).resolve(strict=False)),
        "method_a": a_summary.get("method", args.name_a),
        "method_b": b_summary.get("method", args.name_b),
        "num_queries": n,
        "agreement": agreement,
        "disagreement": float(1.0 - agreement) if n else float("nan"),
        "both_visible": both_visible,
        "both_blocked": both_blocked,
        "a_visible_b_blocked": a_only,
        "a_blocked_b_visible": b_only,
        "a_visible_ratio": float(a_labels.mean()) if n else float("nan"),
        "b_visible_ratio": float(b_labels.mean()) if n else float("nan"),
        "a_query_time_s": a_time,
        "b_query_time_s": b_time,
        "speedup_b_over_a": float(a_time / b_time) if b_time > 0 and math.isfinite(a_time) else float("nan"),
        "a_throughput_qps": a_summary.get("throughput_qps", float("nan")),
        "b_throughput_qps": b_summary.get("throughput_qps", float("nan")),
        "a_peak_rss_gb": a_summary.get("peak_rss_gb", float("nan")),
        "b_peak_rss_gb": b_summary.get("peak_rss_gb", float("nan")),
    }
    print(json.dumps(row, ensure_ascii=False, indent=2))
    if args.output:
        append_csv_row(Path(args.output), row)
        print(f"Comparison row appended to {args.output}")
    return 0


def compare_labels(args: argparse.Namespace) -> int:
    a = load_label_array(Path(args.a))
    b = load_label_array(Path(args.b))
    if a.shape != b.shape:
        raise ValueError(f"Label shapes differ: {a.shape} vs {b.shape}")

    n = int(a.size)
    both_visible = int(np.count_nonzero(a & b))
    both_blocked = int(np.count_nonzero(~a & ~b))
    a_only = int(np.count_nonzero(a & ~b))
    b_only = int(np.count_nonzero(~a & b))
    agree = both_visible + both_blocked
    row = {
        "label_a": args.name_a,
        "label_b": args.name_b,
        "file_a": str(Path(args.a).resolve(strict=False)),
        "file_b": str(Path(args.b).resolve(strict=False)),
        "num_labels": n,
        "agreement": float(agree / n) if n else float("nan"),
        "disagreement": float((a_only + b_only) / n) if n else float("nan"),
        "both_visible": both_visible,
        "both_blocked": both_blocked,
        "a_visible_b_blocked": a_only,
        "a_blocked_b_visible": b_only,
        "a_visible_ratio": float(a.mean()) if n else float("nan"),
        "b_visible_ratio": float(b.mean()) if n else float("nan"),
    }
    print(json.dumps(row, ensure_ascii=False, indent=2))
    if args.output:
        append_csv_row(Path(args.output), row)
        print(f"Agreement row appended to {args.output}")
    return 0


def print_bbox(args: argparse.Namespace) -> int:
    obj_paths = discover_obj_paths(args.obj, args.obj_dir, args.obj_glob)
    bbox_min, bbox_max, vertex_count, elapsed = scan_obj_bbox(obj_paths)
    row = {
        "num_obj_files": len(obj_paths),
        "vertex_count": vertex_count,
        "bbox_min_x": float(bbox_min[0]),
        "bbox_min_y": float(bbox_min[1]),
        "bbox_min_z": float(bbox_min[2]),
        "bbox_max_x": float(bbox_max[0]),
        "bbox_max_y": float(bbox_max[1]),
        "bbox_max_z": float(bbox_max[2]),
        "scan_time_s": elapsed,
    }
    print(json.dumps(row, ensure_ascii=False, indent=2))
    if args.output:
        append_csv_row(Path(args.output), row)
        print(f"BBox row appended to {args.output}")
    return 0


def add_obj_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--obj", action="append", default=[], help="OBJ file, directory, or glob. May be repeated.")
    parser.add_argument("--obj-dir", help="Directory containing OBJ tiles.")
    parser.add_argument("--obj-glob", default="*.obj", help="OBJ search pattern inside directories.")


def add_run_args(parser: argparse.ArgumentParser) -> None:
    add_obj_args(parser)
    parser.add_argument(
        "--targets",
        action="append",
        required=True,
        help="Target points file: .npy, .npz, .csv, or .txt. May be repeated.",
    )
    parser.add_argument("--tx", nargs=3, type=float, required=True, metavar=("X", "Y", "Z"))
    parser.add_argument("--eps", type=float, default=0.01, help="Endpoint tolerance in model units.")
    parser.add_argument("--chunk-size", type=int, default=500_000, help="Number of rays per batch.")
    parser.add_argument("--build-threads", type=int, default=0, help="Open3D scene build threads; 0 means automatic.")
    parser.add_argument("--query-threads", type=int, default=0, help="Open3D cast_rays threads; 0 means automatic.")
    parser.add_argument("--repeat", type=int, default=1, help="Repeat query timing for each target set.")
    parser.add_argument("--output", default="out/bvh_rt_summary.csv", help="Output summary CSV.")
    parser.add_argument("--mesh-report", default="out/bvh_rt_mesh_report.csv", help="Per-OBJ mesh report CSV.")
    parser.add_argument("--labels-dir", help="Optional directory for per-target visibility labels as .npy.")
    parser.add_argument("--results-dir", help="Optional directory for complete visibility result .npz files.")
    parser.add_argument("--result-suffix", default="bvh_rt", help="Suffix for files written under --results-dir.")
    parser.add_argument("--result-include-points", action="store_true", help="Store point coordinates inside result .npz.")
    parser.add_argument("--offset", nargs=3, type=float, metavar=("X", "Y", "Z"), help="Manual offset to subtract.")
    parser.add_argument("--zero-offset", action="store_true", help="Do not shift coordinates.")
    parser.add_argument("--no-mmap", action="store_true", help="Load .npy targets fully instead of memory mapping.")
    parser.add_argument("--progress-every", type=int, default=0, help="Print progress after this many rays; 0 disables.")
    parser.set_defaults(func=run_bvh_rt)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run BVH-RT visibility queries or prepare target/label files."
    )
    subparsers = parser.add_subparsers(dest="command")

    run_parser = subparsers.add_parser("run", help="Run Open3D BVH-RT visibility benchmark.")
    add_run_args(run_parser)

    gen = subparsers.add_parser("generate-targets", help="Generate reproducible local x,y,z target .npy files.")
    add_obj_args(gen)
    gen.add_argument("--bbox", nargs=4, type=float, metavar=("XMIN", "XMAX", "YMIN", "YMAX"))
    gen.add_argument("--bbox-padding", type=float, default=0.0)
    gen.add_argument("--out-dir", default="targets")
    gen.add_argument("--prefix", default="targets")
    gen.add_argument("--spacing", type=float, default=1.0)
    gen.add_argument("--height-range", nargs=3, type=float, default=(10.0, 210.0, 5.0), metavar=("START", "STOP", "STEP"))
    gen.add_argument("--height-list", help="Comma-separated explicit heights. Overrides --height-range.")
    gen.add_argument("--sizes", default=DEFAULT_SIZES)
    gen.add_argument("--dtype", choices=["float32", "float64"], default="float64")
    gen.add_argument("--seed", type=int, default=20260517)
    gen.add_argument("--sample-smallest", action=argparse.BooleanOptionalAction, default=True)
    gen.add_argument("--chunk-rows", type=int, default=1_000_000)
    gen.add_argument("--overwrite", action="store_true")
    gen.add_argument("--manifest", default="out/bvh_rt_targets_manifest.csv")
    gen.set_defaults(func=generate_targets)

    cand = subparsers.add_parser("generate-candidates", help="Generate one fixed-height local candidate grid file.")
    cand.add_argument("--bbox", nargs=4, type=float, required=True, metavar=("XMIN", "XMAX", "YMIN", "YMAX"))
    cand.add_argument("--z", type=float, required=True, help="Fixed slice height in local coordinates.")
    cand.add_argument("--spacing", type=float, default=1.0)
    cand.add_argument("--output", required=True)
    cand.add_argument("--dtype", choices=["float32", "float64"], default="float64")
    cand.add_argument("--manifest", default="out/candidate_manifest.csv")
    cand.set_defaults(func=generate_candidates)

    cmp_parser = subparsers.add_parser("compare-labels", help="Compare two boolean visibility label arrays.")
    cmp_parser.add_argument("--a", required=True, help="First .npy/.npz/.csv label file.")
    cmp_parser.add_argument("--b", required=True, help="Second .npy/.npz/.csv label file.")
    cmp_parser.add_argument("--name-a", default="BVH-RT")
    cmp_parser.add_argument("--name-b", default="VRPF")
    cmp_parser.add_argument("--output", default="out/bvh_rt_label_agreement.csv")
    cmp_parser.set_defaults(func=compare_labels)

    merge_parser = subparsers.add_parser("merge-labels", help="Write labels plus optional performance into a result .npz.")
    merge_parser.add_argument("--points", required=True, help="Common candidate point file.")
    merge_parser.add_argument("--labels", required=True, help="Visibility labels from BVH-RT or VRPF.")
    merge_parser.add_argument("--method", required=True, help="Method name, e.g. BVH-RT or VRPF.")
    merge_parser.add_argument("--output", required=True, help="Output .npz result file.")
    merge_parser.add_argument("--query-time-s", type=float, default=float("nan"))
    merge_parser.add_argument("--peak-rss-gb", type=float, default=float("nan"))
    merge_parser.add_argument("--threads", type=int, default=0)
    merge_parser.add_argument("--summary-csv", default="out/visibility_result_summary.csv")
    merge_parser.add_argument("--include-points", action="store_true")
    merge_parser.add_argument("--no-mmap", action="store_true")
    merge_parser.set_defaults(func=merge_visibility_result)

    conv_parser = subparsers.add_parser("local-to-lonlat", help="Convert local x,y,z candidate points to lon,lat,z.")
    conv_parser.add_argument("--points", required=True)
    conv_parser.add_argument("--output", required=True)
    conv_parser.add_argument("--mode", choices=["linear", "pyproj"], default="linear")
    conv_parser.add_argument("--origin-x", type=float, default=0.0)
    conv_parser.add_argument("--origin-y", type=float, default=0.0)
    conv_parser.add_argument("--origin-lon", type=float, default=0.0)
    conv_parser.add_argument("--origin-lat", type=float, default=0.0)
    conv_parser.add_argument("--meters-to-lon", type=float, default=1.0 / 102500.0)
    conv_parser.add_argument("--meters-to-lat", type=float, default=1.0 / 111320.0)
    conv_parser.add_argument("--source-crs", default="EPSG:2326")
    conv_parser.add_argument("--target-crs", default="EPSG:4326")
    conv_parser.add_argument("--projected-offset-x", type=float, default=0.0)
    conv_parser.add_argument("--projected-offset-y", type=float, default=0.0)
    conv_parser.add_argument("--no-mmap", action="store_true")
    conv_parser.set_defaults(func=convert_local_to_lonlat)

    viz_parser = subparsers.add_parser("visualize-slice", help="Visualize one fixed-height visibility slice as PNG.")
    viz_parser.add_argument("--points", required=True, help="Common candidate point file.")
    viz_parser.add_argument("--result", required=True, help="Result .npz or label file.")
    viz_parser.add_argument("--output", required=True)
    viz_parser.add_argument("--z", type=float)
    viz_parser.add_argument("--z-tol", type=float, default=1e-6)
    viz_parser.add_argument("--spacing", type=float)
    viz_parser.add_argument("--visible-color", nargs=3, type=int, default=(253, 231, 37), metavar=("R", "G", "B"))
    viz_parser.add_argument("--blocked-color", nargs=3, type=int, default=(68, 1, 84), metavar=("R", "G", "B"))
    viz_parser.add_argument("--no-mmap", action="store_true")
    viz_parser.set_defaults(func=visualize_slice)

    result_cmp_parser = subparsers.add_parser("compare-results", help="Compare two method result .npz files.")
    result_cmp_parser.add_argument("--a", required=True)
    result_cmp_parser.add_argument("--b", required=True)
    result_cmp_parser.add_argument("--name-a", default="BVH-RT")
    result_cmp_parser.add_argument("--name-b", default="VRPF")
    result_cmp_parser.add_argument("--output", default="out/visibility_result_comparison.csv")
    result_cmp_parser.set_defaults(func=compare_results)

    bbox_parser = subparsers.add_parser("bbox", help="Scan OBJ vertex bounding box.")
    add_obj_args(bbox_parser)
    bbox_parser.add_argument("--output", help="Optional CSV output path.")
    bbox_parser.set_defaults(func=print_bbox)

    return parser


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = build_parser()
    known_commands = {
        "run",
        "generate-targets",
        "generate-candidates",
        "compare-labels",
        "merge-labels",
        "local-to-lonlat",
        "visualize-slice",
        "compare-results",
        "bbox",
        "-h",
        "--help",
    }
    if argv and argv[0] not in known_commands:
        argv = ["run", *argv]
    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help()
        raise SystemExit(2)
    return args


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
