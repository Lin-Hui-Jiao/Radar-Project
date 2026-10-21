# Fixed-Height 1km Experiment Run Guide

This document explains how to build and run the fixed-height visibility comparison program for the 1km dataset.

## 1. Files Used

The program target is:

```bash
build/bin/fixed_height_experiment
```

The default 1km index file is:

```bash
test_area_1_1.3idx
```

The default output/cache directory is:

```bash
fixed_height_1km_outputs
```

Important: the program default index path is `../../test_area_1_1.3idx`, so the default run command assumes your current working directory is `build/bin`.

## 2. Build

Run from the repository root:

```bash
cd /home/rootjiao/radar-demo-0927-shiyan
cmake --build build --target fixed_height_experiment -j 2
```

If the `build` directory has not been configured yet:

```bash
cd /home/rootjiao/radar-demo-0927-shiyan
cmake -S . -B build
cmake --build build --target fixed_height_experiment -j 2
```


## 3. Recommended Run Method

Run from `build/bin` so the default index path resolves correctly:

```bash
cd /home/rootjiao/radar-demo-0927-shiyan/build/bin
./fixed_height_experiment
```

This uses the defaults:

| Parameter | Default |
| --- | --- |
| Index file | `../../test_area_1_1.3idx` |
| Output directory | `fixed_height_1km_outputs` |
| Radar lon/lat | center of the 1km dataset bounds |
| Radar relative height | `80m` above terrain |
| Fixed-height plane | `50m` ASL |
| Grid resolution | `1m` |
| DEM resolutions | `1m`, `5m`, `10m`, `15m`, `20m` |

## 4. Run With Explicit Radar Settings

Example:

```bash
cd /home/rootjiao/radar-demo-0927-shiyan/build/bin
./fixed_height_experiment \
  --index ../../test_area_1_1.3idx \
  --output-dir fixed_height_1km_outputs \
  --radar-lon 114.1621448335 \
  --radar-lat 22.2783642355 \
  --radar-alt 40 \
  --fixed-height 15
```

`--radar-alt` is relative height above local terrain. The program queries the terrain height from the RTree and converts it to absolute altitude internally.

`--fixed-height` is absolute height above sea level.

## 5. Run From Repository Root

If you run the executable from the repository root, do not rely on the default index path. Pass paths explicitly:

```bash
cd /home/rootjiao/radar-demo-0927-shiyan
./build/bin/fixed_height_experiment \
  --index test_area_1_1.3idx \
  --output-dir build/bin/fixed_height_1km_outputs \
  --radar-alt 80 \
  --fixed-height 50
```

## 6. Useful Options

Show all options:

```bash
./fixed_height_experiment --help
```

Common options:

| Option | Meaning |
| --- | --- |
| `--index <path>` | RTree index file path |
| `--output-dir <dir>` | DEM cache and experiment output directory |
| `--radar-lon <deg>` | Radar longitude |
| `--radar-lat <deg>` | Radar latitude |
| `--radar-alt <m>` | Radar height above terrain |
| `--base-phi <deg>` | Radar base pitch angle |
| `--base-theta <deg>` | Radar base azimuth angle |
| `--fixed-height <m>` | Fixed-height plane, absolute height ASL |
| `--min-lon`, `--min-lat`, `--max-lon`, `--max-lat` | Override dataset bounds |
| `--index-range-x`, `--index-range-y`, `--min-x`, `--min-y` | Override RTree coordinate conversion parameters |

## 7. Output Files

By default, outputs are written under:

```bash
/home/rootjiao/radar-demo-0927-shiyan/build/bin/fixed_height_1km_outputs
```

Generated files:

```bash
1km_dem_1m.tif
1km_dem_5m.tif
1km_dem_10m.tif
1km_dem_15m.tif
1km_dem_20m.tif

1km_fixed_height_result_1m.csv
1km_fixed_height_result_5m.csv
1km_fixed_height_result_10m.csv
1km_fixed_height_result_15m.csv
1km_fixed_height_result_20m.csv

1km_fixed_height_viz_1m.png
1km_fixed_height_viz_5m.png
1km_fixed_height_viz_10m.png
1km_fixed_height_viz_15m.png
1km_fixed_height_viz_20m.png
```

## 8. Cache Warning

The program reuses existing DEM files if the same file name already exists in the output directory.

If you change any of these parameters, use a new `--output-dir` or delete the old DEM files first:

```bash
--min-lon
--min-lat
--max-lon
--max-lat
--index-range-x
--index-range-y
--min-x
--min-y
```

For example:

```bash
cd /home/rootjiao/radar-demo-0927-shiyan/build/bin
rm -f fixed_height_1km_outputs/1km_dem_*.tif
```

## 9. Thread Count

The mesh/DEM comparison loop uses OpenMP. You can control thread count with `OMP_NUM_THREADS`:

```bash
cd /home/rootjiao/radar-demo-0927-shiyan/build/bin
OMP_NUM_THREADS=8 ./fixed_height_experiment
```

## 10. Windows PowerShell Example

From Windows PowerShell:

```powershell
wsl -d Ubuntu-22.04 -- bash -lc "cd /home/rootjiao/radar-demo-0927-shiyan/build/bin && ./fixed_height_experiment"
```

With explicit options:

```powershell
wsl -d Ubuntu-22.04 -- bash -lc "cd /home/rootjiao/radar-demo-0927-shiyan/build/bin && ./fixed_height_experiment --index ../../test_area_1_1.3idx --output-dir fixed_height_1km_outputs --radar-alt 80 --fixed-height 50"
```
