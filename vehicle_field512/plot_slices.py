#!/usr/bin/env python3
"""Slice figures of a vehicle_field512 field, for checking it by eye.

    python3 vehicle_field512/plot_slices.py V1

Inputs:  data/v512_<id>_vd.{f32,json}; data/top_surface_512.f32 (site_check heights, same window); the true mesh
         cross-sections of the slice planes (site_check section, cached as data/section_<plane>.txt).
Outputs: figures/<id>_horizontal.png  horizontal slices with the mesh cross-section of the slice plane
         figures/<id>_vertical.png    vertical slices through the radar (east-west, north-south) with the mesh
                                      cross-section, plus close-ups around the radar
         The top surface (Getheight3d, highest triangle of a column) hides overhangs and recesses, so the figures draw
         the intersection of every triangle with the slice plane instead. Node (i, j, k) is drawn centred on its
         coordinates (E0 + i, N0 + j, altitude k).
         figures/<id>_physics.png     S * R^2 against elevation angle (main lobe constant, side lobes below 0 deg) and
                                      S against R, recomputed from the kernel formula (coordinates by the system PROJ,
                                      cs2cs, i.e. the same library as field_export). This uses the kernel's own
                                      local frame (sphere of radius 6371193 m), so it shows that the formula was
                                      evaluated as written, not that R is the true distance: see phys_check.py
Prints a few numeric checks (non-zero nodes under the top surface, agreement with the formula).
"""
import json
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
import subprocess

D = os.path.dirname(os.path.abspath(__file__))
SID = sys.argv[1] if len(sys.argv) > 1 else "V1"
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#d4d3cf"
RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]   # sequential blue, low -> high
OCC = "#a3a19a"   # occluded (value 0), neutral
VMIN, VMAX = 1e-6, 1e1

for f in ("/mnt/c/Windows/Fonts/msyh.ttc", "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"):
    if os.path.exists(f):
        font_manager.fontManager.addfont(f)
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name()
        break
plt.rcParams.update({"axes.unicode_minus": False, "axes.edgecolor": GRID, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.labelcolor": INK2, "font.size": 9})

pre = os.path.join(D, "data", f"v512_{SID}_vd")
J = json.load(open(pre + ".json"))
g, rd = J["grid"], J["provenance"]["radar"]
nx, ny, nz = g["nx"], g["ny"], g["nz"]
E0, N0, Z0 = g["origin"]
F = np.memmap(pre + ".f32", dtype="<f4", mode="r", shape=(nz, ny, nx))
H = np.fromfile(os.path.join(D, "data", "top_surface_512.f32"), np.float32).reshape(ny, nx).astype(float)
cov = H > -1e8   # columns with triangles (Getheight3d misses the 125 m sub-tile seam lines; see README section 7)
ri, rj, rz = rd["epsg2326"][0] - E0, rd["epsg2326"][1] - N0, rd["abs_alt_m"] - Z0
cov_i = np.nonzero(cov.any(axis=0))[0]; cov_j = np.nonzero(cov.any(axis=1))[0]
cmap = LinearSegmentedColormap.from_list("seq_blue", RAMP)
cmap.set_bad(OCC)
norm = LogNorm(VMIN, VMAX, clip=True)
figdir = os.path.join(D, "figures"); os.makedirs(figdir, exist_ok=True)


def show(ax, img, extent):
    a = np.ma.masked_where(img == 0, np.clip(img, VMIN, VMAX))
    return ax.imshow(a, origin="lower", cmap=cmap, norm=norm, extent=extent, interpolation="nearest", aspect="equal")


IDX = os.path.join(D, "..", "test_area.3idx")
OFF = J["provenance"]["rtree_offset"]
SC = os.path.join(D, "build", "site_check")


def mesh_section(plane, c, a0, a1, b0, b1):
    """Segments (a0, b0, a1, b1) where the triangles cut the plane (site_check section), cached in data/."""
    f = os.path.join(D, "data", f"section_{plane}_{c:g}_{a0:g}_{a1:g}_{b0:g}_{b1:g}.txt")
    if not os.path.exists(f):
        subprocess.run([SC, "section", IDX, repr(OFF[0]), repr(OFF[1]), plane, repr(c), repr(a0), repr(a1), repr(b0),
                        repr(b1), f], check=True, capture_output=True)
    if os.path.getsize(f) == 0:   # plane above every triangle
        return np.zeros((0, 4))
    return np.loadtxt(f, ndmin=2).reshape(-1, 4)


def draw_section(ax, seg, da, db, lw):
    ax.add_collection(LineCollection([((s[0] - da, s[1] - db), (s[2] - da, s[3] - db)) for s in seg], colors=INK,
                                     linewidths=lw, zorder=4))


def radar_mark(ax, x, y):
    ax.plot(x, y, "o", ms=8, mfc="white", mec=INK, mew=1.8, zorder=5)


def legend_items():
    return [Patch(facecolor=OCC, edgecolor="none", label="0（被遮挡）"),
            Line2D([], [], color=INK, lw=0.8, label="网格截面（三角面与切面的交线）"),
            Line2D([], [], color=INK2, lw=1, ls="--", label="索引覆盖范围（框外 6 m 外扩带无几何）"),
            Line2D([], [], marker="o", ls="none", mfc="white", mec=INK, mew=1.8, ms=8, label="雷达")]


title = (f"{SID}：雷达 E {rd['epsg2326'][0]:.0f} / N {rd['epsg2326'][1]:.0f}（窗口内 ({ri:.0f}, {rj:.0f})），"
         f"天线海拔 {rd['abs_alt_m']:.2f} m，俯仰 {rd['elevation_range_deg'][0]:g}°～{rd['elevation_range_deg'][1]:g}°")

# ---------------------------------------------------------------- horizontal slices
levels = [6, int(math.floor(rz)), 30, 60, 120, 250]   # second panel: the layer just below the antenna
fig, axs = plt.subplots(2, 3, figsize=(15.5, 10.4), dpi=130, constrained_layout=True)
fig.patch.set_facecolor(SURFACE)
for ax, z in zip(axs.flat, levels):
    sl = np.asarray(F[z - int(Z0)])
    im = show(ax, sl, (-0.5, nx - 0.5, -0.5, ny - 0.5))
    draw_section(ax, mesh_section("h", Z0 + z, E0, E0 + nx - 1, N0, N0 + ny - 1), E0, N0, 0.5)
    ax.set_xlim(-0.5, nx - 0.5); ax.set_ylim(-0.5, ny - 0.5)
    ax.add_patch(Rectangle((cov_i[0] - 0.5, cov_j[0] - 0.5), cov_i[-1] - cov_i[0] + 1, cov_j[-1] - cov_j[0] + 1, fill=False,
                           ec=INK2, lw=1, ls="--"))
    radar_mark(ax, ri, rj)
    vis = (sl != 0).mean() * 100
    note = "（低于雷达水平面：旁瓣区）" if z < rz else ""
    ax.set_title(f"海拔 {z} m{note}　可见 {vis:.1f}%", color=INK, fontsize=10, loc="left")
    ax.set_xlabel(f"东向 i (m)，E = {E0:.0f} + i"); ax.set_ylabel(f"北向 j (m)，N = {N0:.0f} + j")
cb = fig.colorbar(im, ax=axs, shrink=0.8, pad=0.01)
cb.set_label("功率密度 S（W/m²，对数）", color=INK2); cb.outline.set_visible(False)
fig.legend(handles=legend_items(), loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.45, -0.035))
fig.suptitle("水平切面　" + title, color=INK, fontsize=12, x=0.01, ha="left")
fig.savefig(os.path.join(figdir, f"{SID}_horizontal.png"), facecolor=SURFACE, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- vertical slices through the radar
I0, J0 = int(round(ri)), int(round(rj))
ew = np.asarray(F[:, J0, :])     # [z][i]
ns = np.asarray(F[:, :, I0])     # [z][j]
fig, axs = plt.subplots(2, 2, figsize=(14, 12.5), dpi=130, constrained_layout=True,
                        gridspec_kw={"height_ratios": [1.0, 0.62]})
fig.patch.set_facecolor(SURFACE)
zoom_top = max(60, 10 * math.ceil((rz + 45) / 10))
for col, (sec, seg, a0, ra, name, xl) in enumerate(
        ((ew, mesh_section("ew", N0 + J0, E0, E0 + nx - 1, Z0 - 1, Z0 + nz), E0, ri,
          f"东西向竖直切面（N = {N0 + J0:.0f}，j = {J0}）", "东向 i (m)"),
         (ns, mesh_section("ns", E0 + I0, N0, N0 + ny - 1, Z0 - 1, Z0 + nz), N0, rj,
          f"南北向竖直切面（E = {E0 + I0:.0f}，i = {I0}）", "北向 j (m)"))):
    for row, (x0, x1, z1) in enumerate(((-0.5, sec.shape[1] - 0.5, Z0 + nz - 0.5), (ra - 60, ra + 60, zoom_top))):
        ax = axs[row, col]
        im = show(ax, sec, (-0.5, sec.shape[1] - 0.5, Z0 - 0.5, Z0 + nz - 0.5))
        draw_section(ax, seg, a0, 0.0, 0.6 if row == 0 else 0.9)
        radar_mark(ax, ra, rz)
        ax.axhline(rz, color=INK2, lw=0.8, ls=":")
        ax.set_xlim(max(-0.5, x0), min(sec.shape[1] - 0.5, x1)); ax.set_ylim(Z0 - 0.5, z1)
        ax.set_xlabel(xl); ax.set_ylabel("海拔 (m)")
        ax.set_title(name if row == 0 else name + "　雷达附近放大（点线 = 雷达水平面）", color=INK, fontsize=10, loc="left")
cb = fig.colorbar(im, ax=axs, shrink=0.7, pad=0.01)
cb.set_label("功率密度 S（W/m²，对数）", color=INK2); cb.outline.set_visible(False)
fig.legend(handles=[legend_items()[0], legend_items()[1], legend_items()[3]],
           loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.45, -0.025))
fig.suptitle("竖直切面　" + title, color=INK, fontsize=12, x=0.01, ha="left")
fig.savefig(os.path.join(figdir, f"{SID}_vertical.png"), facecolor=SURFACE, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- physics check (system PROJ + kernel formula)
PI, EM = math.pi, 6371193.0
K = 1.0 * 10 ** (25 / 10) / (4 * PI)          # Pt [W] * G / (4 pi) = 25.165
dl, NEL = 0.006 / (300.0 / 24000.0), 256
def afac(off):   # |sin(N u) / (N sin u)|, u = pi d/lambda sin(off)
    off = np.asarray(off, float); u = PI * dl * np.sin(np.radians(off)); den = NEL * np.sin(u)
    with np.errstate(divide="ignore", invalid="ignore"):
        v = np.abs(np.sin(NEL * u) / den)
    return np.where((off <= 0) | (np.abs(den) < 1e-300), 1.0, v)
rlon, rlat = rd["lonlat"]; ralt = rd["abs_alt_kernel_m"]
rng = np.random.default_rng(11)
ii = rng.integers(0, nx, 400000); jj = rng.integers(0, ny, 400000); kk = rng.integers(0, nz, 400000)
kk[:150000] = rng.integers(0, 40, 150000)   # oversample the band near the radar
kb = int(math.floor(rz))                     # every non-zero node of the layers below the radar horizon
lowz, lowj, lowi = np.nonzero(np.asarray(F[:kb + 1]))
ii = np.concatenate([ii, lowi]); jj = np.concatenate([jj, lowj]); kk = np.concatenate([kk, lowz])
S = np.asarray(F[kk, jj, ii], float)
m = S > 0
ii, jj, kk, S = ii[m], jj[m], kk[m], S[m]
# EPSG:2326 -> EPSG:4326 with the system PROJ (authority axis order: northing easting -> lat lon). pyproj bundles another
# PROJ version whose results differ by ~4e-10 deg (0.04 mm); harmless in the main lobe, but the side lobes are so
# sensitive to the angle that it shows there as up to ~100 % near the nulls.
_out = subprocess.run(["cs2cs", "-f", "%.15f", "EPSG:2326", "EPSG:4326"], check=True, capture_output=True, text=True,
                      input="".join(f"{N0 + j:.6f} {E0 + i:.6f}\n" for i, j in zip(ii, jj))).stdout.split()
lat = np.array(_out[0::3], float); lon = np.array(_out[1::3], float)
rl1, rl2 = math.radians(rlat), np.radians(lat)
x = EM * (rl2 - rl1); y = EM * 2 * np.arcsin(math.cos(rl1) * np.sin(np.radians(rlon - lon) / 2)); z = Z0 + kk - ralt
R2 = x * x + y * y + z * z; R = np.sqrt(R2); phi = np.degrees(np.arcsin(z / R))
pred = (K * afac(-phi) / R2).astype(np.float32).astype(float)   # the file stores float32
rel = np.abs(S - pred) / pred
up = phi >= 0
print(f"formula consistency (same local frame as the kernel; see phys_check.py for true distances): {S.size} non-zero sampled nodes; max |S - formula| / formula: main lobe ({up.sum()}) "
      f"{rel[up].max():.2e} (R > 2 m: {rel[up & (R > 2)].max():.2e}), below the horizon ({(~up).sum()}) {rel[~up].max():.2e} "
      f"; identical after float32 rounding {np.mean(S == pred) * 100:.2f} %")
fig, axs = plt.subplots(1, 2, figsize=(14, 5.2), dpi=130, constrained_layout=True)
fig.patch.set_facecolor(SURFACE)
ax = axs[0]
pp = np.linspace(-90, 90, 20001)
ax.plot(pp, K * afac(-pp), color=INK2, lw=0.35, alpha=0.7, zorder=1, label="公式：25.165 × 方向图因子（俯仰 < 0° 为旁瓣）")
ax.scatter(phi, S * R2, s=7, color=RAMP[3], edgecolors=RAMP[6], linewidths=0.3, zorder=3,
           label=f"文件中的非零节点（主瓣抽样 {up.sum()} 个，水平面以下全部 {(~up).sum()} 个）")
ax.set_yscale("log"); ax.set_xlim(-90, 90); ax.set_ylim(1e-6, 1e2)
ax.axvline(0, color=INK2, lw=0.8, ls=":")
ax.set_xlabel("俯仰角 (°)，相对雷达"); ax.set_ylabel("S × R²（W）")
ax.set_title(f"S × R² 随俯仰角：0°～90° 为常数 25.165；0° 以下为旁瓣（与公式最大相对偏差 {rel[~up].max():.1e}，见说明）", color=INK, fontsize=10, loc="left")
ax.legend(frameon=False, loc="lower right")
ax.grid(True, color=GRID, lw=0.5)
ax = axs[1]
up = phi >= 0
ax.plot(R[up], S[up], "o", ms=2.5, mfc=RAMP[3], mec="none", zorder=3, label="俯仰 ≥ 0° 的非零节点（抽样）")
rr = np.logspace(np.log10(max(0.3, R[up].min())), np.log10(R.max()), 200)
ax.plot(rr, K / rr ** 2, color=INK, lw=0.8, zorder=4, label="25.165 / R²")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("到雷达的距离 R (m)"); ax.set_ylabel("S（W/m²）")
ax.set_title(f"主瓣区 S 与 R（R 为生成器的球面近似距离）：与公式最大相对偏差 {rel[up & (R > 2)].max():.1e}（R > 2 m），{rel[up].max():.1e}（含雷达 2 m 内）",
             color=INK, fontsize=10, loc="left")
ax.legend(frameon=False, loc="upper right")
ax.grid(True, color=GRID, lw=0.5)
fig.suptitle("公式一致性核对（与生成器同一坐标近似；按真实距离的核对见 phys_check.py）　" + title, color=INK, fontsize=12, x=0.01, ha="left")
fig.savefig(os.path.join(figdir, f"{SID}_physics.png"), facecolor=SURFACE, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- numeric checks on the geometry
under = 0; under_cols = 0
k_idx = np.arange(nz)[:, None]
for j0 in range(0, ny, 64):
    blk = np.asarray(F[:, j0:j0 + 64, :])                     # [z][j][i]
    h = H[j0:j0 + 64, :]
    below = (Z0 + np.arange(nz))[:, None, None] < (h[None, :, :] - 0.5)   # clearly under the top surface
    below &= cov[j0:j0 + 64, :][None, :, :]
    nzb = (blk != 0) & below
    under += int(nzb.sum()); under_cols += int(nzb.any(axis=0).sum())
print(f"geometry check: non-zero nodes more than 0.5 m under the top surface (covered columns): {under} "
      f"in {under_cols} columns (possible only under overhangs and in recesses: bridges, canopies, eaves, undercut "
      f"facades; the top surface is the highest triangle of a column, compare with the mesh cross-sections)")
for name in ("horizontal", "vertical", "physics"):
    print("wrote", os.path.join(figdir, f"{SID}_{name}.png"))
