# radar-demo-0927-shiyan 工程说明（每次会话自动读入）

## 交流
用简体中文回复。

## 本工程是什么
硕士论文**第一工作点**（VRPF）：基于倾斜摄影城市网格（R 树 + 线段三角形求交）计算雷达直射功率密度。

## 第二工作点在别处
压缩存储与查询的工作已迁到独立工程 **`/home/rootjiao/pshc`**，背景、约束、方案、进度都在它的 `CLAUDE.md` 和 `docs/` 中。做第二工作点相关的事，先读 `/home/rootjiao/pshc/CLAUDE.md`。

本工程在第二工作点中**只负责生成数据**：
- 待新增 `field_export/` 目录，输出交接格式（`<名称>.f32` + `<名称>.json`，规格见 pshc 的 `CLAUDE.md`）；
- 可以参考 `rise_review/gen_mesh3d.cpp`、`rise_review/gen_sample1024c.cpp`、`rise_review/gen_sample1024d.cpp`。

## 规则
- **不得修改本工程已有的文件**，只能新增文件或目录。工作分支为 `feature/storage-compression`。
- `.gitignore` 会忽略所有 `*.py` 和 `*.md`，提交这类新文件时需要 `git add -f`。

## 生成数据时必须注意
- 雷达绝对高度 = R 树在雷达处的地面高程 + 安装高度（500 m 区域为 164.454871 m，不是打印出来的 164.2 m）。
- CityGuardRadar 的俯仰扫描范围在 `include/hiradar/radar_pool.hpp` 中是 −90°～0°；第二工作点假设全向扫描，导出时用 −90°～90°（做成参数，不改原文件）。
- 经纬度以单精度保存（`DecimalToDMS`、`DMSToDecimal`），会带来约 0.16% 的锯齿噪声；导出工具应同时支持 V-float（与现有程序一致）和 V-double。
- `Radar::Move` 的坐标常量（−1975 / 43）与 1 km、2 km 数据集配置（`index_range_x/y`）换算后相差约 0.36 m / 0.43 m，导出时统一。
- `build/bin/dem_20m.tif` 实际分辨率约 0.52 m。
