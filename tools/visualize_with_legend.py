#!/usr/bin/env python3
"""
雷达功率密度可视化脚本（带完整图例）
使用与C++ fixed_height_comparator.cpp完全相同的Viridis配色和对数映射

特性：
- X轴：经度范围
- Y轴：纬度范围
- Colorbar：功率密度（W/m²），对数刻度
- 与C++代码的颜色映射完全一致
- 保持原始网格精细度（如500×500）

用法：
    python visualize_with_legend.py <csv_file> [radar_lon] [radar_lat]

示例：
    python visualize_with_legend.py fixed_height_result_5m.csv
    python visualize_with_legend.py fixed_height_result_5m.csv 114.1670 22.2806
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import cm
import sys
import os

# Viridis配色方案 - 与C++代码完全一致
VIRIDIS_COLORMAP = [
    (0.00,  68,   1,  84),  # 深紫色
    (0.10,  71,  15, 107),  # 深紫蓝
    (0.20,  72,  36, 117),  # 紫蓝色
    (0.30,  67,  55, 125),  # 蓝紫色
    (0.40,  56,  75, 126),  # 深蓝色
    (0.50,  42,  95, 125),  # 蓝色
    (0.60,  32, 115, 117),  # 蓝绿色
    (0.70,  34, 135,  97),  # 绿蓝色
    (0.80,  68, 154,  68),  # 黄绿色
    (0.90, 122, 170,  40),  # 浅黄绿
    (1.00, 253, 231,  37),  # 亮黄色
]

def create_viridis_colormap():
    """创建Viridis colormap（与C++一致）"""
    colors = []
    for val, r, g, b in VIRIDIS_COLORMAP:
        colors.append([r/255.0, g/255.0, b/255.0])

    cmap = mcolors.LinearSegmentedColormap.from_list('viridis_custom', colors, N=256)
    return cmap

def map_density_to_color_log(density, min_density, max_density):
    """
    功率密度到颜色映射（与C++代码完全一致）

    对数映射 + Gamma校正 (gamma=0.6)
    """
    if density <= 0:
        return np.array([128, 128, 128], dtype=np.uint8)  # 灰色

    if max_density <= min_density or max_density <= 0 or min_density <= 0:
        return np.array([128, 128, 128], dtype=np.uint8)

    # 对数归一化
    log_density = np.log(density)
    log_min = np.log(min_density)
    log_max = np.log(max_density)

    if abs(log_max - log_min) < 1e-6:
        # 返回中间值颜色
        mid_idx = len(VIRIDIS_COLORMAP) // 2
        return np.array(VIRIDIS_COLORMAP[mid_idx][1:], dtype=np.uint8)

    # 归一化到[0, 1]
    normalized = (log_density - log_min) / (log_max - log_min)
    normalized = np.clip(normalized, 0.0, 1.0)

    # Gamma校正 (gamma = 0.6)
    gamma = 0.6
    normalized = np.power(normalized, gamma)

    # 在colormap中查找并插值
    for i in range(len(VIRIDIS_COLORMAP) - 1):
        v1, r1, g1, b1 = VIRIDIS_COLORMAP[i]
        v2, r2, g2, b2 = VIRIDIS_COLORMAP[i + 1]

        if v1 <= normalized <= v2:
            # 线性插值
            t = (normalized - v1) / (v2 - v1)
            r = int(r1 + t * (r2 - r1))
            g = int(g1 + t * (g2 - g1))
            b = int(b1 + t * (b2 - b1))
            return np.array([r, g, b], dtype=np.uint8)

    # 返回最大值颜色
    return np.array(VIRIDIS_COLORMAP[-1][1:], dtype=np.uint8)

def load_and_process_data(csv_file):
    """加载CSV数据并处理"""
    print(f"Loading: {csv_file}")
    df = pd.read_csv(csv_file)

    # 获取网格信息
    unique_lons = sorted(df['lon'].unique())
    unique_lats = sorted(df['lat'].unique())
    n_lon = len(unique_lons)
    n_lat = len(unique_lats)

    print(f"  Grid size: {n_lon} × {n_lat} points")
    print(f"  Lon range: [{unique_lons[0]:.6f}, {unique_lons[-1]:.6f}]")
    print(f"  Lat range: [{unique_lats[0]:.6f}, {unique_lats[-1]:.6f}]")

    # 创建2D网格
    power_grid = np.zeros((n_lat, n_lon))
    visible_grid = np.zeros((n_lat, n_lon), dtype=bool)

    for _, row in df.iterrows():
        i_lon = unique_lons.index(row['lon'])
        i_lat = unique_lats.index(row['lat'])
        power_grid[i_lat, i_lon] = row['power_density']
        visible_grid[i_lat, i_lon] = row['visible_dem']

    # 计算功率密度范围
    visible_powers = power_grid[visible_grid & (power_grid > 0)]

    if len(visible_powers) == 0:
        print("  Warning: No valid power data!")
        min_power = 1e-6
        max_power = 1e-5
    else:
        min_power = visible_powers.min()
        max_power = visible_powers.max()
        print(f"  Power range: {min_power:.3e} - {max_power:.3e} W/m²")
        print(f"  DEM visible points: {visible_grid.sum()} / {visible_grid.size}")

    return df, unique_lons, unique_lats, power_grid, visible_grid, min_power, max_power

def create_image_with_cpp_colormap(power_grid, visible_grid, min_power, max_power):
    """使用C++相同的颜色映射创建图像"""
    n_lat, n_lon = power_grid.shape
    image = np.zeros((n_lat, n_lon, 3), dtype=np.uint8)

    print("  Applying color mapping (C++ algorithm)...")

    for i_lat in range(n_lat):
        for i_lon in range(n_lon):
            if visible_grid[i_lat, i_lon]:
                density = power_grid[i_lat, i_lon]
                color = map_density_to_color_log(density, min_power, max_power)
                # 垂直翻转以匹配图像坐标
                image[n_lat - 1 - i_lat, i_lon] = color
            else:
                # DEM不可见：灰色
                image[n_lat - 1 - i_lat, i_lon] = [128, 128, 128]

    return image

def plot_with_full_legend(image, extent, min_power, max_power, output_file, radar_pos=None):
    """生成带完整图例的图像"""
    fig, ax = plt.subplots(figsize=(14, 10))

    # 显示图像
    im = ax.imshow(image, extent=extent, origin='lower', aspect='auto', interpolation='nearest')

    # 创建自定义colormap用于colorbar
    cmap = create_viridis_colormap()

    # 对数归一化
    norm = mcolors.LogNorm(vmin=min_power, vmax=max_power)

    # 添加colorbar
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04, aspect=20)
    cbar.set_label('Power Density (W/m²)', fontsize=14, fontweight='bold',
                   rotation=270, labelpad=25)

    # 设置colorbar刻度为科学计数法
    cbar.ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:.2e}'))

    # 设置坐标轴
    ax.set_xlabel('Longitude (°)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Latitude (°)', fontsize=14, fontweight='bold')
    ax.set_title('Radar Power Density Distribution (DEM-Based Visibility)\nViridis Colormap with Log Scale',
                 fontsize=16, fontweight='bold', pad=20)

    # 格式化坐标轴刻度
    ax.ticklabel_format(useOffset=False, style='plain')
    ax.tick_params(labelsize=11)

    # 标注雷达位置
    if radar_pos is not None:
        ax.plot(radar_pos[0], radar_pos[1], 'r*', markersize=20,
                markeredgecolor='white', markeredgewidth=2.5,
                label='Radar Position', zorder=10)
        ax.legend(loc='upper right', fontsize=12, framealpha=0.95)

    # 添加网格
    ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5, color='gray')

    # 添加说明文本
    stats_text = 'Color Mapping:\n'
    stats_text += f'• Purple → Cyan → Green → Yellow\n'
    stats_text += f'  (Viridis: perceptually uniform)\n'
    stats_text += f'• Range: {min_power:.2e} - {max_power:.2e} W/m²\n'
    stats_text += f'• Gray: DEM occluded (not visible)\n'
    stats_text += f'• Scale: Logarithmic + Gamma(0.6)'

    ax.text(0.02, 0.98, stats_text,
            transform=ax.transAxes, fontsize=10,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.95,
                     edgecolor='black', linewidth=1))

    # 保存高质量图像
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {output_file}")
    plt.close()

def main():
    if len(sys.argv) < 2:
        print("Usage: python visualize_with_legend.py <csv_file> [radar_lon] [radar_lat]")
        print("\nExample:")
        print("  python visualize_with_legend.py fixed_height_result_5m.csv")
        print("  python visualize_with_legend.py fixed_height_result_5m.csv 114.1670 22.2806")
        sys.exit(1)

    csv_file = sys.argv[1]

    if not os.path.exists(csv_file):
        print(f"Error: File not found: {csv_file}")
        sys.exit(1)

    # 雷达位置（可选）
    radar_pos = None
    if len(sys.argv) >= 4:
        radar_pos = (float(sys.argv[2]), float(sys.argv[3]))
        print(f"Radar position: lon={radar_pos[0]}, lat={radar_pos[1]}")

    print("\n" + "="*70)
    print("Radar Power Density Visualization with Full Legend")
    print("Using Viridis colormap (matching C++ implementation)")
    print("="*70 + "\n")

    # 加载和处理数据
    df, unique_lons, unique_lats, power_grid, visible_grid, min_power, max_power = \
        load_and_process_data(csv_file)

    # 创建图像（使用C++相同的颜色映射算法）
    print("\nCreating image with C++-matching color mapping...")
    image = create_image_with_cpp_colormap(power_grid, visible_grid, min_power, max_power)

    # 图像范围
    extent = [unique_lons[0], unique_lons[-1], unique_lats[0], unique_lats[-1]]

    # 生成输出文件名
    base_name = csv_file.replace('.csv', '')
    output_file = f"{base_name}_with_legend.png"

    # 生成带完整图例的图像
    print("\nGenerating visualization with full legend...")
    plot_with_full_legend(image, extent, min_power, max_power, output_file, radar_pos)

    print("\n" + "="*70)
    print("Visualization completed!")
    print("="*70)
    print(f"\nGenerated file:")
    print(f"  {output_file}")
    print(f"\nFeatures:")
    print(f"  ✓ X-axis: Longitude range")
    print(f"  ✓ Y-axis: Latitude range")
    print(f"  ✓ Colorbar: Power density (W/m²) with log scale")
    print(f"  ✓ Viridis colormap (perceptually uniform, colorblind-friendly)")
    print(f"  ✓ Gray areas: DEM occluded regions")
    print(f"  ✓ High resolution: 300 DPI (publication quality)")
    print(f"\nColor mapping algorithm:")
    print(f"  - Log normalization: handles multiple orders of magnitude")
    print(f"  - Gamma correction (0.6): enhances mid-range contrast")
    print(f"  - Matches C++ implementation exactly\n")

if __name__ == "__main__":
    main()
