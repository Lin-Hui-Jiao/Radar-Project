# 多 OBJ 场景下 BVH-RT 三维遮挡检测实验操作说明

本文档进一步细化审稿补充实验中 **BVH-based 3D ray casting baseline** 的实际操作流程。当前方案只采用一种实现路线：

```text
Python + Open3D RaycastingScene + 多 OBJ 三角网格 + 3D ray casting
```

目标是完成一个真正 3D 的基准方法：

1. 输入多个 `.obj` 文件；
2. 所有 `.obj` 文件共同组成同一个局部坐标系下的三维城市模型；
3. 输入一个雷达点 `P_tx` 和大量目标点 `P_rx`;
4. 对每条 `P_tx -> P_rx` 三维线段进行遮挡检测；
5. 输出目标点是否可见、核心查询耗时、吞吐量、可见比例和内存；
6. 用同一批目标点与 VRPF 进行公平对比。

已经新增脚本：

```text
tools/bvh_rt_multiobj_experiment.py
```

该脚本就是本实验的建议实现骨架。

## 1. 方法定义

论文中可将该基准方法命名为：

```text
BVH-RT
```

完整含义：

```text
BVH-based three-dimensional ray casting baseline
```

它不是 DSM/DTM，不是 2.5D 方法，而是真正使用 OBJ 三角网格的 3D 方法。

核心判定是：

```text
给定雷达点 P_tx 和目标点 P_rx：

1. 构造方向向量 d = normalize(P_rx - P_tx)
2. 从 P_tx 向 P_rx 发射三维 ray
3. 查询 ray 是否在到达 P_rx 之前与任意 OBJ 三角面相交
4. 如果相交，则 blocked = true
5. 如果不相交，则 visible = true
```

由于 ray 理论上是无限长的，所以必须限制交点距离：

```text
blocked = isfinite(t_hit) and t_hit < distance(P_tx, P_rx) - 2 * eps
```

其中 `eps` 是端点容差，用来避免发射点或目标点附近的自相交。

## 2. 多 OBJ 的处理方式

### 2.1 推荐方式

对于由多个 OBJ 文件组成的城市模型，推荐做法是：

```text
把所有 OBJ tile 逐个读入，然后加入同一个 Open3D RaycastingScene。
```

也就是说：

```text
OBJ tile 1 \
OBJ tile 2  \
OBJ tile 3   ---> one RaycastingScene ---> batch ray casting
...          /
OBJ tile N  /
```

不要逐个 OBJ 文件分别查询再合并结果。原因是：

1. 每条 ray 只需要查询一次全场景；
2. 查询逻辑更接近“同一个三维城市模型”；
3. 避免多个 OBJ 查询结果合并时的计时和逻辑偏差；
4. 更适合在论文中叙述为 conventional 3D mesh-based baseline。

### 2.2 坐标要求

所有 OBJ 文件必须满足：

```text
同一个局部坐标系
同一个单位
同一个 z-up 高度方向
```

雷达点和目标点也必须在同一个坐标系下。

如果 OBJ 文件是局部坐标系，例如：

```text
x, y 表示局部平面坐标
z 表示高度
单位为米
```

那么输入点也必须是：

```text
P_tx = (x, y, z)
P_rx = (x, y, z)
```

不能直接输入经纬度。

## 3. 安装运行环境

建议使用 Python 3.9 或更高版本。

安装依赖：

```powershell
pip install open3d numpy psutil
```

其中：

- `open3d` 用于读取 OBJ 和构建 RaycastingScene；
- `numpy` 用于批量目标点和 ray 数组；
- `psutil` 用于统计进程内存；如果不安装，也可以运行，但内存字段会为空。

## 4. 输入文件组织

建议目录结构如下：

```text
project/
  obj_tiles/
    tile_001.obj
    tile_002.obj
    tile_003.obj
    ...
  targets/
    targets_1e5.npy
    targets_5e5.npy
    targets_1e6.npy
    targets_5e6.npy
    targets_1e7.npy
  tools/
    bvh_rt_multiobj_experiment.py
  out/
```

其中 `obj_tiles/` 存放所有 OBJ 文件。

## 5. 目标点文件格式

脚本支持以下目标点格式：

```text
.npy
.npz
.csv
.txt
```

推荐使用 `.npy`，因为 1000 万点 CSV 会非常大，读取也慢。

目标点数组形状必须为：

```text
N x 3
```

每一行是一个目标点：

```text
x, y, z
```

CSV 示例：

```csv
x,y,z
100.0,200.0,30.0
101.0,200.0,30.0
102.0,200.0,30.0
```

`.npy` 示例：

```python
import numpy as np

targets = np.array([
    [100.0, 200.0, 30.0],
    [101.0, 200.0, 30.0],
    [102.0, 200.0, 30.0],
], dtype=np.float64)

np.save("targets_1e5.npy", targets)
```

## 6. 建议的目标点生成方案

用户计划测试：

```text
100,000
500,000
1,000,000
5,000,000
10,000,000
```

建议采用“固定水平高度切片”的方式生成目标点，便于论文解释，也便于和现有 VRPF 网格实验对应。

### 6.1 基本设计

假设研究区是：

```text
500 m x 500 m
```

平面采样间隔为：

```text
1 m
```

则每个高度切片有：

```text
500 x 500 = 250,000 points
```

不同规模可通过抽样或叠加高度切片生成：

| Query set | Number of target points | Suggested construction |
|---|---:|---|
| Q1 | 100,000 | 从一个高度切片中固定抽样 100,000 点 |
| Q2 | 500,000 | 2 个高度切片 |
| Q3 | 1,000,000 | 4 个高度切片 |
| Q4 | 5,000,000 | 20 个高度切片 |
| Q5 | 10,000,000 | 40 个高度切片 |

高度可设为：

```text
10 m, 15 m, 20 m, ..., 205 m
```

也可以根据论文低空 UAV 场景设置为：

```text
20 m, 25 m, 30 m, ..., 215 m
```

最重要的是：**BVH-RT 和 VRPF 必须使用完全相同的目标点文件。**

### 6.2 为什么不推荐随机点作为主实验

随机点也可以，但需要额外说明：

1. 随机种子；
2. 采样范围；
3. 如何排除建筑物内部点；
4. 如何保证不同方法点集一致；
5. 随机点是否覆盖了代表性高度层。

固定高度切片更适合审稿补实验，因为逻辑简单、可复现、和原文网格实验衔接自然。

## 7. 生成目标点的示例脚本

可以新建一个临时脚本生成 `.npy` 文件。下面只是示例，需要把范围改成你的 OBJ 局部坐标范围。

```python
import numpy as np
from pathlib import Path

out_dir = Path("targets")
out_dir.mkdir(exist_ok=True)

# 修改为你的研究区局部坐标范围
x_min, x_max = 0.0, 500.0
y_min, y_max = 0.0, 500.0
spacing = 1.0

xs = np.arange(x_min, x_max, spacing)
ys = np.arange(y_min, y_max, spacing)
xx, yy = np.meshgrid(xs, ys, indexing="xy")
xy = np.column_stack([xx.ravel(), yy.ravel()])

heights = np.arange(10.0, 210.0, 5.0)  # 40 slices

all_points = []
for z in heights:
    zz = np.full((xy.shape[0], 1), z)
    pts = np.column_stack([xy, zz])
    all_points.append(pts)

points = np.vstack(all_points).astype(np.float64)

# 固定抽样 100,000 点
rng = np.random.default_rng(20260517)
idx_1e5 = rng.choice(points.shape[0], size=100_000, replace=False)
np.save(out_dir / "targets_1e5.npy", points[idx_1e5])

np.save(out_dir / "targets_5e5.npy", points[:500_000])
np.save(out_dir / "targets_1e6.npy", points[:1_000_000])
np.save(out_dir / "targets_5e6.npy", points[:5_000_000])
np.save(out_dir / "targets_1e7.npy", points[:10_000_000])
```

如果研究区不是 `0-500`，必须改成实际局部坐标范围，例如：

```text
x_min = OBJ bbox min x
x_max = OBJ bbox max x
y_min = OBJ bbox min y
y_max = OBJ bbox max y
```

## 8. 脚本运行命令

### 8.1 单个规模测试

假设 OBJ 文件夹为：

```text
D:\paper\paper\paper1\obj_tiles
```

目标点文件为：

```text
D:\paper\paper\paper1\targets\targets_1e5.npy
```

雷达点为：

```text
P_tx = (250.0, 250.0, 30.0)
```

运行：

```powershell
python tools\bvh_rt_multiobj_experiment.py `
  --obj-dir D:\paper\paper\paper1\obj_tiles `
  --targets D:\paper\paper\paper1\targets\targets_1e5.npy `
  --tx 250.0 250.0 30.0 `
  --eps 0.01 `
  --chunk-size 500000 `
  --query-threads 32 `
  --output out\bvh_rt_summary.csv `
  --mesh-report out\bvh_rt_mesh_report.csv
```

### 8.2 一次跑五个规模

```powershell
python tools\bvh_rt_multiobj_experiment.py `
  --obj-dir D:\paper\paper\paper1\obj_tiles `
  --targets D:\paper\paper\paper1\targets\targets_1e5.npy `
  --targets D:\paper\paper\paper1\targets\targets_5e5.npy `
  --targets D:\paper\paper\paper1\targets\targets_1e6.npy `
  --targets D:\paper\paper\paper1\targets\targets_5e6.npy `
  --targets D:\paper\paper\paper1\targets\targets_1e7.npy `
  --tx 250.0 250.0 30.0 `
  --eps 0.01 `
  --chunk-size 500000 `
  --query-threads 32 `
  --repeat 3 `
  --output out\bvh_rt_summary.csv `
  --mesh-report out\bvh_rt_mesh_report.csv `
  --labels-dir out\bvh_rt_labels
```

说明：

- `--repeat 3` 表示每个规模重复 3 次；
- `--labels-dir` 会保存每个目标点的可见/遮挡标签，后续可与 VRPF 输出做一致率比较；
- 如果 1000 万点太耗时，可以先设置 `--repeat 1`。

### 8.3 多个 OBJ 路径的其他写法

可以直接指定多个 OBJ：

```powershell
python tools\bvh_rt_multiobj_experiment.py `
  --obj tile_001.obj `
  --obj tile_002.obj `
  --obj tile_003.obj `
  --targets targets\targets_1e6.npy `
  --tx 250 250 30
```

也可以指定通配符：

```powershell
python tools\bvh_rt_multiobj_experiment.py `
  --obj "obj_tiles\*.obj" `
  --targets targets\targets_1e6.npy `
  --tx 250 250 30
```

但最推荐的是：

```powershell
--obj-dir obj_tiles
```

## 9. 输出文件说明

### 9.1 summary 文件

默认输出：

```text
out/bvh_rt_summary.csv
```

主要字段：

| Field | Meaning |
|---|---|
| `method` | 方法名，固定为 `BVH-RT-Open3D` |
| `target_file` | 当前目标点文件 |
| `repeat` | 重复实验编号 |
| `num_obj_files` | OBJ 文件数量 |
| `total_vertices` | 总顶点数 |
| `total_triangles` | 总三角面数 |
| `tx_x, tx_y, tx_z` | 雷达点 |
| `eps` | 端点容差 |
| `chunk_size` | 每批 ray 数量 |
| `mesh_read_time_s` | OBJ 读取时间 |
| `mesh_clean_shift_time_s` | mesh 清理和平移时间 |
| `scene_build_time_s` | RaycastingScene 构建时间 |
| `num_queries` | 查询数量 |
| `query_time_s` | 核心遮挡查询时间 |
| `throughput_qps` | 每秒查询条数 |
| `num_visible` | 可见目标数 |
| `num_blocked` | 遮挡目标数 |
| `visible_ratio` | 可见比例 |
| `peak_rss_gb` | 查询阶段峰值内存 |
| `labels_file` | 可见/遮挡标签文件 |

论文中核心使用：

```text
num_queries
query_time_s
throughput_qps
visible_ratio
peak_rss_gb
```

索引构建成本单独报告：

```text
mesh_read_time_s
scene_build_time_s
```

### 9.2 mesh report 文件

默认输出：

```text
out/bvh_rt_mesh_report.csv
```

它记录每个 OBJ tile 的：

```text
obj_file
vertices
triangles
read_time_s
clean_shift_time_s
add_to_scene_time_s
```

这可以用于论文补充说明：

```text
The multi-tile OBJ scene contains X OBJ files and Y triangular facets.
```

### 9.3 labels 文件

若使用：

```text
--labels-dir out\bvh_rt_labels
```

脚本会输出：

```text
targets_1e6_rep1_visible.npy
```

每个 `.npy` 是一个布尔数组：

```text
True  = visible
False = blocked
```

后续可与 VRPF 的输出标签比较：

```python
agreement = np.mean(bvh_visible == vrpf_visible)
```

## 10. 与 VRPF 对比时必须统一的内容

为了公平，BVH-RT 和 VRPF 必须统一：

| Item | Requirement |
|---|---|
| OBJ 数据 | 同一批 OBJ 文件 |
| 雷达点 | 完全相同的 `P_tx` |
| 目标点 | 完全相同的 `.npy` 或 CSV |
| 坐标系 | 同一个局部坐标系 |
| 单位 | 同一单位，建议米 |
| 遮挡定义 | 两点之间是否存在三角面相交 |
| 端点容差 | 尽量统一 `eps` |
| 计时范围 | 核心遮挡查询时间单独报告 |
| 硬件 | 同一台机器 |
| 线程数 | 若比较并行版本，使用同样线程数 |

不要把 BVH-RT 的 OBJ 读取时间和 VRPF 的纯查询时间直接比较。建议在论文中拆开：

```text
1. preprocessing/index construction time
2. repeated visibility-query time
```

## 11. 推荐实验矩阵

建议最终实验表如下：

| Method | Threads | 100k | 500k | 1M | 5M | 10M |
|---|---:|---:|---:|---:|---:|---:|
| BVH-RT | 32 | time | time | time | time | time |
| VRPF-1T | 1 | time | time | time | time | time |
| VRPF-Full | 32 | time | time | time | time | time |

另外记录：

| Method | Index build time | Peak memory | Notes |
|---|---:|---:|---|
| BVH-RT | TBD | TBD | multi-OBJ Open3D RaycastingScene |
| VRPF | TBD | TBD | proposed R-tree index |

若篇幅有限，论文正文只放核心查询时间表，把索引构建和内存放补充材料或简表。

## 12. 论文叙述模板

### 12.1 方法段

```latex
\revised{To provide a truly three-dimensional baseline, we implemented a conventional BVH-based ray-casting method, denoted as BVH-RT. The baseline uses the same multi-tile OBJ model as VRPF. All OBJ tiles are loaded in the same local Cartesian coordinate frame and inserted into one ray-casting scene. For each radar-target pair, a 3D ray is cast from the radar position to the target point. If the first ray-mesh intersection occurs before the target point, the path is classified as obstructed; otherwise, it is classified as visible.}
```

### 12.2 多 OBJ 说明

```latex
\revised{The study area is represented by multiple OBJ files. In the BVH-RT baseline, these OBJ tiles are not queried independently. Instead, all triangular facets from the tiles are inserted into a single ray-casting scene, so that each transmitter-target segment is tested once against the complete 3D city model.}
```

### 12.3 计时说明

```latex
\revised{For fairness, the reported query runtime includes only the repeated line-of-sight obstruction tests. OBJ loading, mesh cleaning, and BVH/R-tree construction are reported separately as reusable preprocessing costs. The same transmitter position, target points, local coordinate system, and hardware platform are used for BVH-RT and VRPF.}
```

### 12.4 规模说明

```latex
\revised{Five query scales were evaluated: \(10^5\), \(5\times10^5\), \(10^6\), \(5\times10^6\), and \(10^7\) radar-target segments. The target points were generated from horizontal slices of the same study area. The exact same target-point files were used by BVH-RT and VRPF.}
```

### 12.5 结果解释

根据实际结果填入：

```latex
\revised{The BVH-RT baseline and VRPF use the same 3D geometry and therefore produce highly consistent obstruction labels. The comparison shows that the advantage of VRPF is not merely due to using a 3D model. Under the same 3D input and target-point sets, VRPF reduces the repeated visibility-query time relative to the conventional BVH-RT baseline, while the single-threaded VRPF variant further helps separate the contribution of the proposed query structure from that of OpenMP parallelization.}
```

如果 VRPF-1T 没有明显快于 BVH-RT，要改成更稳妥的表述：

```latex
\revised{The single-threaded VRPF variant achieves comparable query efficiency to BVH-RT, while the full VRPF implementation improves throughput through point-wise parallelization.}
```

## 13. 审稿回复模板

```text
Thank you for this constructive suggestion. In the revised manuscript, we added a truly three-dimensional baseline, denoted as BVH-RT. This baseline uses the same multi-tile OBJ scene model as VRPF. All OBJ files are loaded in the same local Cartesian coordinate frame and inserted into one ray-casting scene. For each transmitter-target pair, BVH-RT casts a 3D ray and determines obstruction by checking whether the ray intersects any triangular facet before reaching the target point.

We evaluated \(10^5\), \(5\times10^5\), \(10^6\), \(5\times10^6\), and \(10^7\) radar-target segments using the same radar position and target-point files as VRPF. OBJ loading and BVH/R-tree construction are reported separately, while the main runtime comparison focuses on the repeated visibility-query stage. We also included a small ablation study with BVH-RT, single-threaded VRPF, and the full OpenMP-parallel VRPF implementation. These additions clarify that the benefit of VRPF is not solely attributable to the use of 3D geometry, but also to the proposed query organization and parallel execution strategy.
```

## 14. 实验执行顺序建议

建议按以下顺序推进，避免一开始就跑 1000 万点浪费时间：

1. 准备 OBJ 文件夹，确认所有 OBJ 处于同一局部坐标系。
2. 选定雷达点 `P_tx`，确认它在场景合理位置。
3. 生成一个很小的测试目标点文件，例如 1000 点。
4. 运行 BVH-RT 脚本，确认能成功加载多个 OBJ 并输出 visible/blocked。
5. 运行 100,000 点，检查耗时和可见比例是否合理。
6. 保存 labels，与 VRPF 对同一批点的结果计算一致率。
7. 运行 500,000 和 1,000,000 点。
8. 若内存和耗时可接受，再运行 5,000,000 和 10,000,000 点。
9. 对 100k、500k、1M 尽量重复 3 次；对 5M、10M 可根据时间重复 1 到 3 次。
10. 整理 `bvh_rt_summary.csv`，填入论文表格。

## 15. 常见问题

### 15.1 运行很慢怎么办

优先调整：

```text
--chunk-size 200000
--chunk-size 500000
--chunk-size 1000000
```

chunk 太小会增加调度开销，chunk 太大会增加内存压力。

### 15.2 内存不足怎么办

降低：

```text
--chunk-size
```

例如：

```powershell
--chunk-size 100000
```

也可以暂时不保存 labels：

```text
不要设置 --labels-dir
```

### 15.3 坐标值太大怎么办

脚本默认会扫描所有 OBJ 的包围盒中心，并对 mesh 和输入点同时减去同一个 offset。这样不改变几何关系，但能提高 `float32` ray casting 的数值稳定性。

除非你很确定坐标值很小，否则不要使用：

```text
--zero-offset
```

### 15.4 为什么会有少量 BVH-RT 与 VRPF 不一致

常见原因：

1. 端点容差不同；
2. 目标点靠近墙面或地面；
3. 退化三角形处理不同；
4. Open3D 使用 float32 ray，VRPF 使用 double；
5. 两边是否排除了建筑内部点。

如果不一致率很低，例如低于 0.5%，可以在论文中说明差异主要位于几何边界附近。若不一致率很高，则需要检查坐标系、单位和 `eps`。

## 16. 最终应交付的实验结果

至少应得到：

1. `out/bvh_rt_summary.csv`;
2. `out/bvh_rt_mesh_report.csv`;
3. BVH-RT 在五个规模下的查询时间；
4. BVH-RT 的吞吐量；
5. BVH-RT 的可见比例；
6. BVH-RT 的峰值内存；
7. BVH-RT 与 VRPF 的标签一致率；
8. BVH-RT、VRPF-1T、VRPF-Full 的消融对比表。

这些结果足以回应审稿人提出的“缺少真正相关 3D 方法对比”的问题。

