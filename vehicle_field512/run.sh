#!/bin/sh
# Vehicle radar fields, 512 x 512 x 512 (see README.md). For every site of sites.tsv (or the given ids):
#   1. site_check site   — the site must be street-level ground with open sky (exit on failure)
#   2. site_check cover  — columns of the window without triangles (recorded in the log)
#   3. field_export      — the field: R-tree visibility (Intersect3d) + amplitude kernel, V-double, elevation 0..90
#   4. field_check stats / random — SHA-256, visibility statistics, 10,000 nodes recomputed independently (bit-exact)
#   5. check_handoff.py  — the files are a valid PSHC input (hand-off format; exit on failure)
#   6. symbolic links data/v512_<id>_vd.{f32,json} in $PSHC_DATA (default /home/rootjiao/pshc/data if it exists;
#      PSHC_DATA= disables), so that pshc_encode can read the field as data/v512_<id>_vd
# Outputs: data/v512_<id>_vd.{f32,json,log}
#   vehicle_field512/run.sh          all rows
#   vehicle_field512/run.sh V1       only V1
set -e
D=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$D/.." && pwd)
IDX=${IDX:-$R/test_area.3idx}
# R-tree local = EPSG:2326 - model_transform translation of the OBJ tiles (config.json of every tile)
OFFX=-836974.635391304153
OFFY=-815456.572217391338
# window: the 500 m index area (E 835000..835500, N 815500..816000) widened by 6 m on every side, altitude 0..511 m
E0=834995; N0=815495; NXY=512; NZ=512; Z0=0
ELEV_MIN=0; ELEV_MAX=90
T=${THREADS:-8}
PSHC_DATA=${PSHC_DATA-/home/rootjiao/pshc/data}
FE=$R/field_export/build
SC=$D/build/site_check
if [ ! -x "$FE/field_export" ] || [ ! -x "$FE/field_check" ]; then
  cmake -S "$R/field_export" -B "$FE" -DCMAKE_BUILD_TYPE=Release > /dev/null && cmake --build "$FE" -j > /dev/null
fi
if [ ! -x "$SC" ]; then
  cmake -S "$D" -B "$D/build" -DCMAKE_BUILD_TYPE=Release > /dev/null && cmake --build "$D/build" -j > /dev/null
fi
mkdir -p "$D/data"
want=" $* "
grep -v '^#' "$D/sites.tsv" | while read -r id E N h note; do
  [ -n "$id" ] || continue
  if [ "$want" != "  " ]; then
    case "$want" in *" $id "*) ;; *) continue ;; esac
  fi
  p=$D/data/v512_${id}_vd
  echo "== $id: E $E N $N, +$h m"
  if ! "$SC" site "$IDX" $OFFX $OFFY "$E" "$N" "$h" 6 > "$p.log" 2>&1; then
    cat "$p.log"; echo "site $id fails the vehicle-site rules; choose another point (see the passing nodes above)"; exit 1
  fi
  "$SC" cover "$IDX" $OFFX $OFFY $E0 $N0 $NXY $NXY >> "$p.log" 2>&1
  "$FE/field_export" --index "$IDX" --grid metric --origin $E0 $N0 --n $NXY $NXY $NZ --spacing 1 1 1 --z0 $Z0 \
    --radar "$E" "$N" --radar-height "$h" --elev-min $ELEV_MIN --elev-max $ELEV_MAX --version v-double \
    --rtree-offset $OFFX $OFFY --out "$p" --threads $T >> "$p.log" 2>&1
  "$FE/field_check" stats "$p" >> "$p.log" 2>&1
  "$FE/field_check" random "$p" "$IDX" 10000 >> "$p.log" 2>&1
  if ! python3 "$D/check_handoff.py" "$p" >> "$p.log" 2>&1; then
    cat "$p.log"; echo "field $id is not a valid PSHC input"; exit 1
  fi
  if [ -n "$PSHC_DATA" ] && [ -d "$PSHC_DATA" ]; then
    ln -sf "$p.f32" "$PSHC_DATA/" && ln -sf "$p.json" "$PSHC_DATA/"
    echo "linked into $PSHC_DATA: $(basename "$p").f32, $(basename "$p").json" >> "$p.log"
  fi
  cat "$p.log"
done
