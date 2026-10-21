# Fixed-Height 1km EPSG:2326 Migration Notes

## Current Decision

The fixed-height experiment now uses a strict EPSG:2326 meter grid instead of an EPSG:4326 degree grid.

| Item | Value |
| --- | --- |
| Index file | `test_area_1_1.3idx` |
| Lon/lat metadata bounds | `min_lon=114.157292795`, `min_lat=22.273848477`, `max_lon=114.166996872`, `max_lat=22.282879994` |
| Projected meter range | `x=[834250, 835250]`, `y=[815000, 816000]` |
| RTree offsets | `index_range_x=-2724.635391`, `index_range_y=-456.572217` |
| Fixed-height grid resolution | `1.0m` |
| DEM resolution set | `1m`, `5m`, `10m`, `15m`, `20m` |

## Implemented Changes

| Area | Change |
| --- | --- |
| Dataset config | Added `area_width_m=1000.0` and `area_height_m=1000.0`. |
| DEM generation | Outputs EPSG:2326 GeoTIFFs with meter pixel size. `1m` DEM is exactly `1000 x 1000` for the default area. |
| Fixed-height grid | Points are generated from EPSG:2326 pixel/grid centers, then converted back to lon/lat for radar power-density calculation. |
| DEM visibility | Added projected-meter visibility checking via `DEMVisibility::IsOccludedProjected()`. |
| Mesh visibility | Uses projected/RTree coordinates directly for segment intersection. |
| CSV output | Exports `lon,lat,x2326,y2326,rtree_x,rtree_y,alt,visible_mesh,visible_dem,disagreement,power_density`. |
| Output naming | Uses `fixed_height_1km_epsg2326_outputs` and `1km_epsg2326_*` file prefixes to avoid WGS84 DEM cache reuse. |

## Not Changed

| Area | Status |
| --- | --- |
| `src/server_new2.cpp` | Not modified. |
| `src/test.cpp` | Not modified. |
| `src/radar.cpp` range logic | Not modified. |
| `include/hiradar/radar_pool.hpp` action distance `300` | Not modified. |
| Radar physical model parameters | Not modified. |

## Validation Commands

```bash
cmake --build build --target fixed_height_experiment -j 2
```

```bash
cd /home/rootjiao/radar-demo-0927-shiyan/build/bin
./fixed_height_experiment --help
```

Small smoke test:

```bash
cd /home/rootjiao/radar-demo-0927-shiyan/build/bin
./fixed_height_experiment \
  --index ../../test_area_1_1.3idx \
  --output-dir fixed_height_epsg2326_smoke_outputs \
  --area-width 2 \
  --area-height 2 \
  --radar-alt 80 \
  --fixed-height 50
```

Full default run:

```bash
cd /home/rootjiao/radar-demo-0927-shiyan/build/bin
./fixed_height_experiment \
  --index ../../test_area_1_1.3idx \
  --output-dir fixed_height_1km_epsg2326_outputs \
  --radar-alt 80 \
  --fixed-height 50
```
