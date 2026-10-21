# field_export：导出单部雷达的三维功率密度场（PSHC 交接格式）

为第二工作点（PSHC 压缩，工程 `/home/rootjiao/pshc`）生成输入数据。本目录是新增的，独立构建，不修改本仓库已有文件。

## 构建

```sh
cmake -S field_export -B field_export/build -DCMAKE_BUILD_TYPE=Release
cmake --build field_export/build -j
```

依赖：g++（C++17、OpenMP）、PROJ。编译选项 `-O2 -ffp-contract=off`（与在 CSV 上逐位核对过的原型相同，禁止 FMA 合并和 fast-math）。

## 计算方法

| 部分 | 实现 |
| --- | --- |
| 遮挡 | `RTree3d::Intersect3d` + `IntersectCallback`（`src/occlusion_utils.cpp`），目标点用双精度坐标 |
| 雷达绝对高度 | 雷达处 `Getheight3d` 的地面高程 + 安装高度（与 `Radar::Move` 相同） |
| 幅值 V-float | `ReplayKernel`：`Radar::CalculateSinglePointPowerDensity` 的 float32 计算链逐行移植，与现有程序逐位一致 |
| 幅值 V-double | `CleanKernel`：同一公式的双精度实现 |
| 俯仰扫描范围 | 参数 `--elev-min/--elev-max`，默认 −90°～90°（全向扫描，无旁瓣），不改 `radar_pool.hpp` |

两个核函数在 `radar_kernels.hpp`，从 `pshc/prototypes/rise_prototype/src/common.hpp` 原样复制。

## 用法

```
field_export --index FILE.3idx --grid lonlat|metric --origin X Y --n NX NY NZ [--spacing SX SY SZ]
             --z0 Z --radar X Y [--radar-height H] [--elev-min DEG] [--elev-max DEG]
             [--version v-float|v-double] [--rtree-offset OX OY] [--lat-ref LAT|LAT0,LAT1]
             --out PREFIX [--threads N]
```

| 参数 | 含义 |
| --- | --- |
| `--grid lonlat` | 经纬度网格（复现 500 m 数据）：`--origin` 为经纬度；`--spacing` 的米数按 111320 m/度和 cos(lat_ref) 换算成度，与 `FixedHeightComparator::RunExperiment` 相同；`--lat-ref A,B` 表示 (A+B)/2，默认取网格中心纬度 |
| `--grid metric` | EPSG:2326 米制网格（正式数据）：`--origin` 为东、北坐标 |
| `--radar` | \|X\| ≤ 360 视为经纬度，否则为 EPSG:2326；V-float 下经纬度雷达先经过 float 度分秒量化（同 `Radar::Move`） |
| `--rtree-offset` | R 树局部坐标 = EPSG:2326 + 偏移。默认 (−836974.6354, −815456.5722)（1 km、2 km 配置）；master `Radar::Move` 为 (−836975, −815457) |

节点 (i, j, k) 的位置：lonlat 为 (lon0 + i·dlon, lat0 + j·dlat, z0 + k·dz)，metric 为 (E0 + i·dx, N0 + j·dy, z0 + k·dz)，高度为绝对高程。

## 输出

- `PREFIX.f32`：小端 float32，无文件头，按 [z][y][x] 排列（x 最快），单位 W/m²，被遮挡为 0。
- `PREFIX.json`：`grid`（nx、ny、nz、dx、dy、dz、crs、origin）、`value`、`sha256`、`provenance`（雷达位置与实际绝对高度、坐标偏移、俯仰范围、数据版本、索引文件、本仓库提交号、命令行、耗时）。
- PSHC 的解码器只用 `grid`；`provenance` 只作记录。

## 校验工具 field_check

```
field_check stats  PREFIX                  # sha256、可见比例、列跳变次数分布
field_check csv    PREFIX CSV [TOL]        # 与 fixed_height_result_*.csv 中两者都非零的点比较（默认容差 1.2e-10）
field_check random PREFIX INDEX [N] [SEED] # 独立重算 N 个随机节点（Intersect3d + 核函数），要求逐位一致
```

## 已验证的运行

```sh
# 500 m 原型数据复现（与 pshc/prototypes/rise_review/gen_mesh3d.cpp 一致：master 坐标偏移，俯仰 −90°～0°）
field_export/build/field_export --index test_area.3idx --grid lonlat --origin 114.164571 22.278365 \
  --n 500 503 160 --lat-ref 22.278365,22.28288 --z0 1 --radar 114.1670 22.2806 --radar-height 80 \
  --version v-float --rtree-offset -836975 -815457 --elev-min -90 --elev-max 0 --out mesh500_vf_e0
field_export/build/field_check stats mesh500_vf_e0
field_export/build/field_check csv mesh500_vf_e0 build/bin/fixed_height_result_1m.csv

# 1024 m 窗口中心低空样本（2 km 索引，雷达在窗口中心、地面以上 80 m，全向扫描）
field_export/build/field_export --index test_area_2_2.3idx --grid metric --origin 836872 814622 \
  --n 256 256 128 --z0 0 --radar 837000 814750 --radar-height 80 --version v-float --out w1024c_vf_256
field_export/build/field_check random w1024c_vf_256 test_area_2_2.3idx 10000
```

结果见 pshc 工程 `docs/实现进度.md` 与阶段 1 总结。
