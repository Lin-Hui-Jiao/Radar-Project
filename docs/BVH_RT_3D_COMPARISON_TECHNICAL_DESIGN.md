# BVH-RT 3D 对比实验技术设计与实现说明

## 1. 背景与目标

本次新增需求是为现有 VRPF 方法补充一个真正三维条件下的对比实验。现有工程已经有两类相关实现：

- `src/radar.cpp` 和 `src/occlusion_utils.cpp`：通过 R-tree 候选三角面筛选，再执行线段-三角形相交判断，作为 VRPF 遮挡检测核心。
- `src/scaling_memory_benchmark.cpp`：已经形成较完整的可扩展性能实验框架，支持随机三维目标点、VRPF/DSM 对比、线程数、核心查询计时和峰值内存记录。

新增的 3D baseline 命名为 `BVH-RT`，实现路线为：

```text
Python + Open3D RaycastingScene + 多 OBJ 三角网格 + 三维 ray casting
```

它不使用 DSM/DTM，也不是 2.5D 栅格方法，而是直接把 OBJ 三角面加入同一个 BVH 场景，然后对 `P_tx -> P_rx` 线段做三维遮挡检测。

## 2. 对现有工程的调研结论

### 2.1 VRPF 当前遮挡语义

现有 C++ 遮挡逻辑的本质是：

1. 将雷达点和目标点转换到 R-tree 使用的局部三维坐标。
2. 用线段端点的三维包围盒检索候选三角面。
3. 对候选三角面执行 Moller-Trumbore 线段相交检测。
4. 一旦命中遮挡三角面，提前终止查询。

关键代码位置：

- `src/radar.cpp`：`Radar::isOccluded`、`CapablePowerDensity`、`CalculateSinglePointPowerDensity`
- `src/occlusion_utils.cpp`：`lineSegIntersectTri`、`IntersectCallback`
- `src/scaling_memory_benchmark.cpp`：`RunVrpf`、`LineSegmentIntersectsTriangleFast`

因此 BVH-RT 的公平对比对象不是能量密度公式，而是同一批线段在三角网格中的 `visible/blocked` 标签和核心遮挡查询时间。

### 2.2 坐标系统要求

当前 VRPF 相关实验里，常见流程是：

```text
WGS84 lon/lat -> EPSG:2326 projected x/y -> subtract min_x/min_y -> add index_range_x/index_range_y
```

而 Open3D 读取 OBJ 后只理解 OBJ 顶点自身坐标。因此 BVH-RT 的输入必须满足：

- OBJ 顶点、雷达点 `--tx`、目标点文件都在同一个局部笛卡尔坐标系。
- 单位一致，建议为米。
- 高度方向一致，默认 `z-up`。

如果目标点来自经纬度，必须先转换成与 OBJ 完全一致的局部坐标后再输入 Python 工具。

## 3. 已实现内容

正式工具位于：

```text
tools/bvh_rt_multiobj_experiment.py
```

根目录保留兼容入口：

```text
bvh_rt_multiobj_experiment.py
```

为了避免 `.gitignore` 忽略新增实验脚本和本文档，已对 `.gitignore` 增加精确例外。

### 3.1 `run`：运行 BVH-RT 查询

功能：

- 支持 `--obj-dir`、重复 `--obj`、OBJ glob。
- 将所有 OBJ tile 加入同一个 `Open3D RaycastingScene`。
- 自动扫描 OBJ bbox 中心并对 mesh、tx、targets 同步平移，降低大坐标转 `float32` ray 时的数值误差。
- 支持 `.npy`、`.npz`、`.csv`、`.txt` 目标点。
- 对 `.npy` 默认使用 memory map，避免 1000 万点一次性读入内存。
- 分 chunk 批量 ray casting。
- 输出 summary CSV、mesh report CSV、可选逐点布尔标签。

核心遮挡判定：

```text
d = normalize(P_rx - P_tx)
ray_origin = P_tx + eps * d
blocked = isfinite(t_hit) and t_hit < distance(P_tx, P_rx) - 2 * eps
visible = not blocked
```

`eps` 同时避开发射点附近和目标点附近的端点自相交。

### 3.2 `generate-targets`：生成固定高度切片目标点

功能：

- 可由手工 bbox 生成：

```text
--bbox XMIN XMAX YMIN YMAX
```

- 也可由 OBJ bbox 自动生成：

```text
--obj-dir obj_tiles
```

- 支持固定高度序列、采样间距、目标规模列表、随机种子。
- 使用 `.npy` memmap 写出，适合 10M 点。

默认规模：

```text
100000,500000,1000000,5000000,10000000
```

默认高度：

```text
10, 15, 20, ..., 205 m
```

### 3.3 `compare-labels`：标签一致率

功能：

- 比较 BVH-RT 与 VRPF 的可见性布尔标签。
- 输出 agreement、disagreement、两边可见比例和四类混淆统计。
- 可追加写入 CSV。

这一步是论文里说明“二者使用相同 3D 几何后标签是否一致”的关键证据。

### 3.4 `bbox`：快速扫描 OBJ bbox

功能：

- 只扫描 OBJ 的 `v x y z` 顶点行。
- 不依赖 Open3D。
- 用于确认研究区局部坐标范围和目标点生成范围。

## 4. 推荐运行流程

### 4.1 安装依赖

建议使用 Python 3.10、3.11 或 3.12：

```bash
python3.10 -m pip install open3d numpy psutil
```

`psutil` 可选；未安装时内存字段为 `nan`。

注意：本机验证时 Windows 侧默认 `python` 为 3.13.5，`pip install open3d` 没有匹配 wheel。因此正式运行 BVH-RT 时应切换到 WSL Python 3.10 或 Conda Python 3.10/3.11/3.12 环境。

### 4.2 检查 OBJ bbox

```bash
python tools/bvh_rt_multiobj_experiment.py bbox \
  --obj-dir obj_tiles
```

根据输出确认 `x/y/z` 范围是否符合本地米制坐标。

### 4.3 生成固定高度目标点

使用 OBJ bbox：

```bash
python tools/bvh_rt_multiobj_experiment.py generate-targets \
  --obj-dir obj_tiles \
  --spacing 1.0 \
  --height-range 10 210 5 \
  --sizes 100000,500000,1000000,5000000,10000000 \
  --out-dir targets \
  --overwrite
```

或使用明确研究区范围：

```bash
python tools/bvh_rt_multiobj_experiment.py generate-targets \
  --bbox 0 500 0 500 \
  --spacing 1.0 \
  --height-range 10 210 5 \
  --sizes 100000,500000,1000000,5000000,10000000 \
  --out-dir targets \
  --overwrite
```

### 4.4 小样例冒烟测试

先准备一个 1000 点目标文件：

```bash
python tools/bvh_rt_multiobj_experiment.py generate-targets \
  --bbox 0 20 0 20 \
  --spacing 1.0 \
  --height-list 10,20,30 \
  --sizes 1000 \
  --out-dir targets_smoke \
  --overwrite
```

再运行 BVH-RT：

```bash
python tools/bvh_rt_multiobj_experiment.py run \
  --obj-dir obj_tiles \
  --targets targets_smoke/targets_1000.npy \
  --tx 250 250 30 \
  --eps 0.01 \
  --chunk-size 100000 \
  --query-threads 8 \
  --labels-dir out/bvh_rt_labels \
  --output out/bvh_rt_summary.csv \
  --mesh-report out/bvh_rt_mesh_report.csv
```

### 4.5 正式五规模实验

```bash
python tools/bvh_rt_multiobj_experiment.py run \
  --obj-dir obj_tiles \
  --targets targets/targets_100000.npy \
  --targets targets/targets_500000.npy \
  --targets targets/targets_1000000.npy \
  --targets targets/targets_5000000.npy \
  --targets targets/targets_10000000.npy \
  --tx 250 250 30 \
  --eps 0.01 \
  --chunk-size 500000 \
  --query-threads 32 \
  --repeat 3 \
  --labels-dir out/bvh_rt_labels \
  --output out/bvh_rt_summary.csv \
  --mesh-report out/bvh_rt_mesh_report.csv
```

### 4.6 与 VRPF 标签比较

假设 VRPF 已输出同一批目标点的布尔可见标签：

```bash
python tools/bvh_rt_multiobj_experiment.py compare-labels \
  --a out/bvh_rt_labels/targets_1000000_rep1_visible.npy \
  --b out/vrpf_labels/targets_1000000_visible.npy \
  --name-a BVH-RT \
  --name-b VRPF \
  --output out/bvh_rt_label_agreement.csv
```

## 5. 输出解释

### 5.1 `out/bvh_rt_summary.csv`

核心字段：

- `num_queries`：目标点数量。
- `query_time_s`：纯 BVH ray casting 查询时间。
- `throughput_qps`：每秒查询条数。
- `num_visible` / `num_blocked`：可见/遮挡数量。
- `visible_ratio`：可见比例。
- `peak_rss_gb`：查询阶段观测到的峰值 RSS。
- `mesh_read_time_s`、`scene_build_time_s`：OBJ 读取和 BVH 场景构建成本，应作为预处理成本单独报告。

### 5.2 `out/bvh_rt_mesh_report.csv`

每个 OBJ tile 的顶点数、三角面数、读取/清理/加入场景时间。论文里可用来说明：

```text
The multi-tile OBJ scene contains X OBJ files and Y triangular facets.
```

### 5.3 标签文件

`--labels-dir` 会写出：

```text
targets_1000000_rep1_visible.npy
```

其中：

```text
True  = visible
False = blocked
```

## 6. 对比实验口径建议

公平对比时必须统一：

- 同一批 OBJ/三角面数据。
- 同一个雷达点。
- 同一批目标点文件。
- 同一局部坐标系和单位。
- 同一硬件平台。
- 线程数分组明确，例如 `BVH-RT-32T`、`VRPF-1T`、`VRPF-Full`。
- 预处理时间和核心查询时间分开报告。

建议论文主表：

| Method | Threads | 100k | 500k | 1M | 5M | 10M |
|---|---:|---:|---:|---:|---:|---:|
| BVH-RT | 32 | query time | query time | query time | query time | query time |
| VRPF-1T | 1 | query time | query time | query time | query time | query time |
| VRPF-Full | 32 | query time | query time | query time | query time | query time |

另附预处理表：

| Method | Build time | Peak memory | Notes |
|---|---:|---:|---|
| BVH-RT | OBJ read + scene build | RSS | Open3D RaycastingScene |
| VRPF | R-tree load/build | RSS | Project R-tree index |

## 7. 主要风险点与修改意见

1. 坐标系统必须先统一。若 OBJ 是局部坐标而 VRPF 输入是经纬度，不能直接比较，需要导出同一局部坐标目标点或给 VRPF 增加同目标点输入模式。
2. 现有 `src/fixed_height_comparator.cpp` 中 mesh 可见性逻辑有被注释的痕迹，不建议作为本次 3D baseline 的参考结果来源；应优先使用 `src/scaling_memory_benchmark.cpp` 的 VRPF 查询路径。
3. Open3D ray 使用 `float32`，VRPF 精确相交主要为 `double`，几何边界附近可能存在少量标签差异。论文中应报告一致率，并说明差异集中于边界/端点容差。
4. 不要把 BVH-RT 的 OBJ 读取和建树时间直接与 VRPF 的纯查询时间比较。两者都应拆分为预处理和重复查询。
5. 10M 点标签文件本身约 10 MB，目标点 `float64` 约 240 MB；建议保留 `.npy`，避免 CSV。
6. 若 OBJ 坐标绝对值很大，不要使用 `--zero-offset`；默认 bbox 中心平移更稳。

## 8. 验证清单

- `python -m py_compile tools/bvh_rt_multiobj_experiment.py`
- `python tools/bvh_rt_multiobj_experiment.py --help`
- `python tools/bvh_rt_multiobj_experiment.py generate-targets --help`
- 小 OBJ + 小目标点运行 `run`，确认 visible/blocked 数量合理。
- 对同一目标点运行 VRPF，使用 `compare-labels` 计算一致率。
