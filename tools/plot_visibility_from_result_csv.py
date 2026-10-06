#!/usr/bin/env python3
"""Plot BVH/VRPF visibility result CSV by lon/lat.

Expected CSV columns include:
point_id, local_x, local_y, local_z, lon, lat, bvh_visible, bvh_blocked
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_visibility_csv(path: Path, visible_column: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lon_values: list[float] = []
    lat_values: list[float] = []
    visible_values: list[bool] = []

    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"CSV has no header: {path}")

        required = {"lon", "lat", visible_column}
        missing = required.difference(reader.fieldnames)
        if missing:
            raise ValueError(f"CSV is missing required columns: {sorted(missing)}")

        for row in reader:
            lon_values.append(float(row["lon"]))
            lat_values.append(float(row["lat"]))
            visible_values.append(int(float(row[visible_column])) != 0)

    if not lon_values:
        raise ValueError(f"CSV contains no data rows: {path}")

    return (
        np.asarray(lon_values, dtype=np.float64),
        np.asarray(lat_values, dtype=np.float64),
        np.asarray(visible_values, dtype=np.bool_),
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot visibility result CSV by longitude/latitude.")
    parser.add_argument("--input", default="out/bvh_rt_candidates_z30_result.csv", help="Input result CSV.")
    parser.add_argument("--output", default="out/bvh_rt_candidates_z30_visibility.png", help="Output PNG path.")
    parser.add_argument("--visible-column", default="bvh_visible", help="Visibility column name.")
    parser.add_argument("--title", default="BVH-RT Visibility at Fixed Height")
    parser.add_argument("--point-size", type=float, default=1.0)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--fig-width", type=float, default=8.0)
    parser.add_argument("--fig-height", type=float, default=7.0)
    parser.add_argument("--radar-lon", type=float, default=None)
    parser.add_argument("--radar-lat", type=float, default=None)
    parser.add_argument("--visible-color", default="#f2d13d")
    parser.add_argument("--blocked-color", default="#2b2f77")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)

    lon, lat, visible = load_visibility_csv(input_path, args.visible_column)
    blocked = ~visible

    num_points = int(visible.size)
    num_visible = int(visible.sum())
    num_blocked = int(blocked.sum())
    visible_ratio = num_visible / num_points if num_points else 0.0

    fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height), constrained_layout=True)

    ax.scatter(
        lon[blocked],
        lat[blocked],
        s=args.point_size,
        c=args.blocked_color,
        marker="s",
        linewidths=0,
        alpha=0.95,
        label=f"Blocked ({num_blocked})",
    )
    ax.scatter(
        lon[visible],
        lat[visible],
        s=args.point_size,
        c=args.visible_color,
        marker="s",
        linewidths=0,
        alpha=0.95,
        label=f"Visible ({num_visible})",
    )

    if args.radar_lon is not None and args.radar_lat is not None:
        ax.scatter(
            [args.radar_lon],
            [args.radar_lat],
            s=55,
            c="#e31a1c",
            marker="*",
            edgecolors="white",
            linewidths=0.7,
            label="Radar",
            zorder=5,
        )

    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(f"{args.title}\nvisible ratio = {visible_ratio:.4f}, points = {num_points}")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, linewidth=0.3, alpha=0.35)
    ax.legend(loc="best", frameon=True)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=args.dpi)
    plt.close(fig)

    print(f"Input: {input_path}")
    print(f"Output: {output_path}")
    print(f"Points: {num_points}")
    print(f"Visible: {num_visible}")
    print(f"Blocked: {num_blocked}")
    print(f"Visible ratio: {visible_ratio:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
