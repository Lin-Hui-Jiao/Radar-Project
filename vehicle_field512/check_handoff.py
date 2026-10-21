#!/usr/bin/env python3
"""Checks that a field is a valid PSHC input (hand-off format, pshc/CLAUDE.md and pshc/src/field.cpp).

    python3 vehicle_field512/check_handoff.py PREFIX [PREFIX ...]

Required by the PSHC reader / encoder (failure -> exit code 1):
  <prefix>.json  grid.nx, ny, nz positive integers; dx, dy, dz > 0; crs; origin of 3 numbers; sha256 (64 hex digits)
  <prefix>.f32   exactly nx * ny * nz little-endian float32, [z][y][x] with x fastest; SHA-256 equal to the .json;
                 every value finite and >= 0 (pshc_encode rejects negative, infinite and NaN values)
Required by the hand-off convention (failure -> exit code 1):
  value: unit W/m^2, dtype float32, endian little, order zyx, occluded 0
  provenance: radar position and absolute altitude, elevation range, data version, R-tree index, generator commit
Reported only: number of non-zero nodes, value range, negative zeros (-0.0 compares equal to 0 in PSHC).
"""
import hashlib
import json
import os
import sys

import numpy as np


def check(prefix):
    errs = []
    js = json.load(open(prefix + ".json"))
    g = js.get("grid", {})
    try:
        nx, ny, nz = (int(g[k]) for k in ("nx", "ny", "nz"))
        if min(nx, ny, nz) <= 0 or any(g[k] != int(g[k]) for k in ("nx", "ny", "nz")):
            errs.append("grid: nx, ny, nz must be positive integers")
    except (KeyError, TypeError, ValueError):
        return [f"grid: nx, ny, nz missing or invalid"]
    for k in ("dx", "dy", "dz"):
        if not (isinstance(g.get(k), (int, float)) and g[k] > 0): errs.append(f"grid: {k} missing or not > 0")
    if not isinstance(g.get("crs"), str): errs.append("grid: crs missing")
    if not (isinstance(g.get("origin"), list) and len(g["origin"]) == 3): errs.append("grid: origin must have 3 numbers")
    v = js.get("value", {})
    for k, want in (("unit", "W/m^2"), ("dtype", "float32"), ("endian", "little"), ("order", "zyx"), ("occluded", 0)):
        if v.get(k) != want: errs.append(f"value.{k} = {v.get(k)!r}, expected {want!r}")
    sha = js.get("sha256", "")
    if not (isinstance(sha, str) and len(sha) == 64 and all(c in "0123456789abcdef" for c in sha)):
        errs.append("sha256 missing or not 64 hex digits")
    pv = js.get("provenance", {}); rd = pv.get("radar", {})
    for k in ("epsg2326", "lonlat", "abs_alt_m", "elevation_range_deg"):
        if k not in rd: errs.append(f"provenance.radar.{k} missing")
    for k in ("data_version", "index", "generator_commit", "rtree_offset"):
        if k not in pv: errs.append(f"provenance.{k} missing")
    n = nx * ny * nz
    size = os.path.getsize(prefix + ".f32")
    if size != 4 * n:
        errs.append(f".f32 has {size} bytes, grid needs {4 * n}")
        return errs
    h = hashlib.sha256()
    with open(prefix + ".f32", "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""): h.update(b)
    if h.hexdigest() != sha: errs.append(f"SHA-256 of the .f32 {h.hexdigest()} differs from the .json")
    F = np.memmap(prefix + ".f32", dtype="<f4", mode="r", shape=(nz, ny, nx))
    bad = negz = nnz = 0
    vmin, vmax = np.inf, 0.0
    for k0 in range(0, nz, 32):
        a = np.asarray(F[k0:k0 + 32])
        bad += int((~np.isfinite(a)).sum() + (a < 0).sum())
        negz += int(((a == 0) & np.signbit(a)).sum())
        p = a[a > 0]
        nnz += p.size
        if p.size: vmin = min(vmin, float(p.min())); vmax = max(vmax, float(p.max()))
    if bad: errs.append(f"{bad} values are negative, infinite or NaN")
    print(f"{prefix}: grid {nx} x {ny} x {nz}, spacing ({g.get('dx')}, {g.get('dy')}, {g.get('dz')}) m, {g.get('crs')}, "
          f"origin {g.get('origin')}; {pv.get('data_version')}, radar alt {rd.get('abs_alt_m')} m, elevation "
          f"{rd.get('elevation_range_deg')} deg")
    print(f"  {size} bytes, sha256 {'ok' if h.hexdigest() == sha else 'MISMATCH'}; non-zero {nnz} ({nnz / n:.4%}), "
          f"S {vmin:.3g}..{vmax:.3g} W/m^2; negative zeros {negz}")
    return errs


rc = 0
for p in sys.argv[1:]:
    e = check(p)
    print("  PSHC input: " + ("OK" if not e else "FAIL"))
    for x in e: print("   - " + x)
    rc |= bool(e)
sys.exit(rc)
