#!/usr/bin/env python3
"""Physics checks of a vehicle_field512 field that do not reuse the generator's geometry.

    python3 vehicle_field512/phys_check.py V1h

plot_slices.py recomputes the kernel formula with the kernel's own local frame (sphere of radius 6371193 m, the
"near the equator" approximation of caltools::lon2xyz), so it only shows that the formula was evaluated as written.
Here the main-lobe amplitude is compared with K / R^2 for the true distance:
  R_grid  Euclidean distance in EPSG:2326 + altitude (the frame of the grid and of the occlusion test)
  R_geod  WGS84 geodesic (pyproj.Geod) horizontally, altitude difference vertically
  R_kern  the kernel's distance (for reference)
and the prediction from the radii of curvature of the WGS84 ellipsoid at the radar latitude is printed. The side-lobe
check compares the stored pattern factor S R^2 / K below the horizon with the array factor |AF| and |AF|^2.
pyproj bundles another PROJ version than the generator (about 0.04 mm apart): irrelevant for these checks, but do not
use it for bit-exact comparisons.
"""
import json
import math
import os
import sys

import numpy as np
from pyproj import Geod, Transformer

D = os.path.dirname(os.path.abspath(__file__))
SID = sys.argv[1] if len(sys.argv) > 1 else "V1"
pre = os.path.join(D, "data", f"v512_{SID}_vd")
J = json.load(open(pre + ".json"))
g, rd = J["grid"], J["provenance"]["radar"]
nx, ny, nz = g["nx"], g["ny"], g["nz"]
E0, N0, Z0 = g["origin"]
F = np.memmap(pre + ".f32", dtype="<f4", mode="r", shape=(nz, ny, nx))
rE, rN = rd["epsg2326"]; rlon, rlat = rd["lonlat"]; ralt = rd["abs_alt_kernel_m"]
K = rd["Pt_kW"] * 1000.0 * 10 ** (rd["G_dB"] / 10) / (4 * math.pi)
EM = 6371193.0                                  # caltools EARTH_RADIUS (km) * 1000
NEL, DL = rd["N"], rd["d_m"] / (300.0 / rd["f_MHz"])
to_ll = Transformer.from_crs("EPSG:2326", "EPSG:4326", always_xy=True)


def kernel_xyz(lon, lat, z):
    x = EM * np.radians(lat - rlat)                                                   # north
    y = EM * 2 * np.arcsin(math.cos(math.radians(rlat)) * np.sin(np.radians(rlon - lon) / 2))   # west
    return x, y, z


def af(off):
    u = math.pi * DL * np.sin(np.radians(off))
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.abs(np.sin(NEL * u) / (NEL * np.sin(u)))


# ---------------------------------------------------------------- main lobe: S R^2 / K with true distances
rng = np.random.default_rng(5)
n = 600000
i = rng.integers(0, nx, n); j = rng.integers(0, ny, n); k = rng.integers(0, nz, n)
S = np.asarray(F[k, j, i], float)
E = E0 + i; N = N0 + j; z = Z0 + k - ralt
m = (S > 0) & (z >= 0)
i, j, k, S, E, N, z = i[m], j[m], k[m], S[m], E[m], N[m], z[m]
lon, lat = to_ll.transform(E, N)
dE, dN = E - rE, N - rN
R_grid = np.sqrt(dE ** 2 + dN ** 2 + z ** 2)
_, _, s = Geod(ellps="WGS84").inv(np.full_like(lon, rlon), np.full_like(lat, rlat), lon, lat)
R_geod = np.sqrt(s ** 2 + z ** 2)
x, y, _ = kernel_xyz(lon, lat, z)
R_kern = np.sqrt(x ** 2 + y ** 2 + z ** 2)
far = R_grid > 20
print(f"{SID}: radar E {rE} N {rN} alt {ralt:.3f} m, K = Pt G / 4 pi = {K:.6f} W; {far.sum()} main-lobe nodes with R > 20 m")
for name, R in (("R_kern (generator frame)", R_kern), ("R_grid (EPSG:2326)", R_grid), ("R_geod (WGS84 geodesic)", R_geod)):
    q = S * R ** 2 / K - 1
    print(f"  S R^2 / K - 1, {name:25s}: {q[far].min():+.4%} .. {q[far].max():+.4%}, mean |.| {np.abs(q[far]).mean():.4%}")
q = R_geod / R_grid - 1
print(f"  R_geod / R_grid - 1: {q[far].min():+.1e} .. {q[far].max():+.1e}")
azg = np.degrees(np.arctan2(dE, dN)) % 360
el = np.degrees(np.arcsin(z / R_grid))
for a0 in range(0, 360, 45):
    sel = (el < 10) & (R_grid > 50) & (np.abs(((azg - a0 + 180) % 360) - 180) < 5)
    if sel.sum():
        print(f"  azimuth {a0:3d} deg, elevation < 10 deg (n = {sel.sum():5d}): median S R_grid^2 / K - 1 = "
              f"{np.median(S[sel] * R_grid[sel] ** 2 / K - 1):+.4%}")
a, e2 = 6378137.0, 6.69437999014e-3
sp2 = math.sin(math.radians(rlat)) ** 2
M = a * (1 - e2) / (1 - e2 * sp2) ** 1.5; Nv = a / math.sqrt(1 - e2 * sp2)
print(f"  WGS84 radii at latitude {rlat:.4f}: meridional {M:.0f} m, prime vertical {Nv:.0f} m (generator sphere {EM:.0f} m)")
print(f"  prediction: north / south {(M / EM) ** 2 - 1:+.4%}, east / west {(Nv / EM) ** 2 - 1:+.4%}, vertical 0; "
      f"anisotropic K = ({K * (Nv / EM) ** 2:.4f}, {K * (M / EM) ** 2:.4f}, {K:.4f}) for (E, N, z)")

# ---------------------------------------------------------------- side lobes: stored pattern factor vs |AF|, |AF|^2
kb = int(math.floor(ralt - Z0))
zz, jj, ii = np.nonzero(np.asarray(F[:kb + 1]))
if len(zz):
    Sb = np.asarray(F[:kb + 1])[zz, jj, ii].astype(float)
    lo, la = to_ll.transform(E0 + ii, N0 + jj)
    x, y, dz = kernel_xyz(lo, la, Z0 + zz - ralt)
    R2 = x * x + y * y + dz * dz
    off = -np.degrees(np.arcsin(dz / np.sqrt(R2)))
    pf = Sb * R2 / K
    oo = np.linspace(0.47, 0.93, 4001); pk = oo[np.argmax(af(oo))]
    sel = (off > 0.6) & (off < 0.85)
    print(f"side lobes: {len(Sb)} non-zero nodes below the horizon; first side-lobe peak at {pk:.4f} deg: "
          f"|AF| = {af(pk):.4f} ({10 * math.log10(af(pk)):.1f} dB as a power factor), |AF|^2 = {af(pk) ** 2:.4f} "
          f"({20 * math.log10(af(pk)):.1f} dB); stored maximum between 0.60 and 0.85 deg: {pf[sel].max() if sel.any() else float('nan'):.4f}")
    print(f"  stored factor / |AF|: median {np.median(pf / af(off)):.6f} (the generator uses |AF|, not |AF|^2)")
