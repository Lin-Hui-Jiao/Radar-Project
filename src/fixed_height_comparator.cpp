#include "hiradar/fixed_height_comparator.hpp"
#include "hiradar/radar.hpp"  // 引入Radar类定义
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <png.h>
#include <proj.h>
#include <chrono>
#include <omp.h>

namespace hiradar {

// 颜色结构
struct Color { unsigned char r, g, b; };

// 颜色表点
struct ColorPoint {
    float value;
    int r, g, b;
};

// Viridis配色方案 - 感知均匀的顺序色谱
// 优点：色盲友好、黑白打印友好、感知均匀、单调递增亮度
// 配色：深紫色 → 深蓝色 → 蓝绿色 → 绿色 → 黄绿色 → 黄色
static const std::vector<ColorPoint> colormap = {
    {0.00f,  68,   1,  84}, // 深紫色 - 最低信号
    {0.10f,  71,  15, 107}, // 深紫蓝
    {0.20f,  72,  36, 117}, // 紫蓝色
    {0.30f,  67,  55, 125}, // 蓝紫色
    {0.40f,  56,  75, 126}, // 深蓝色
    {0.50f,  42,  95, 125}, // 蓝色
    {0.60f,  32, 115, 117}, // 蓝绿色
    {0.70f,  34, 135,  97}, // 绿蓝色
    {0.80f,  68, 154,  68}, // 黄绿色
    {0.90f, 122, 170,  40}, // 浅黄绿
    {1.00f, 253, 231,  37}  // 亮黄色 - 最高信号
};

// 线性插值
static Color lerp(const ColorPoint& c1, const ColorPoint& c2, float value) {
    float t = (value - c1.value) / (c2.value - c1.value);
    unsigned char r = static_cast<unsigned char>(c1.r + t * (c2.r - c1.r));
    unsigned char g = static_cast<unsigned char>(c1.g + t * (c2.g - c1.g));
    unsigned char b = static_cast<unsigned char>(c1.b + t * (c2.b - c1.b));
    return {r, g, b};
}

// 能量密度到颜色映射（对数映射 + Gamma校正）
static Color mapDensityToColorLog(float density, float min_density, float max_density) {
    if (density <= 0) {
        return {128, 128, 128}; // 灰色：无信号
    }

    if (max_density <= min_density || max_density <= 0 || min_density <= 0) {
        return {128, 128, 128};
    }

    // 对数归一化（处理能量跨越多个数量级的情况）
    float log_density = std::log(density);
    float log_min = std::log(min_density);
    float log_max = std::log(max_density);

    if (std::abs(log_max - log_min) < 1e-6) {
        return {static_cast<unsigned char>(colormap[colormap.size()/2].r),
                static_cast<unsigned char>(colormap[colormap.size()/2].g),
                static_cast<unsigned char>(colormap[colormap.size()/2].b)};
    }

    // 归一化到[0, 1]
    float normalized_v = (log_density - log_min) / (log_max - log_min);
    normalized_v = std::max(0.0f, std::min(1.0f, normalized_v));

    // Gamma校正：增强中间值对比度
    // gamma < 1.0 会增强低值区域的对比度
    // gamma > 1.0 会增强高值区域的对比度
    // gamma = 0.6 是一个较好的平衡值
    const float gamma = 0.6f;
    normalized_v = std::pow(normalized_v, gamma);

    // 边界处理
    if (normalized_v <= colormap.front().value) {
        return {static_cast<unsigned char>(colormap.front().r),
                static_cast<unsigned char>(colormap.front().g),
                static_cast<unsigned char>(colormap.front().b)};
    }
    if (normalized_v >= colormap.back().value) {
        return {static_cast<unsigned char>(colormap.back().r),
                static_cast<unsigned char>(colormap.back().g),
                static_cast<unsigned char>(colormap.back().b)};
    }

    // 在colormap中查找并插值
    for (size_t i = 0; i < colormap.size() - 1; ++i) {
        const auto& c1 = colormap[i];
        const auto& c2 = colormap[i + 1];
        if (normalized_v >= c1.value && normalized_v <= c2.value) {
            return lerp(c1, c2, normalized_v);
        }
    }

    return {static_cast<unsigned char>(colormap.back().r),
            static_cast<unsigned char>(colormap.back().g),
            static_cast<unsigned char>(colormap.back().b)};
}
//这是一个工具函数，用于将经度或纬度（度）转换为米，或者将米转换为经度或纬度（度）。
double FixedHeightComparator::DegreesToMeters(double degrees, double latitude, bool is_longitude) {
    const double METERS_PER_DEG_LAT = 111320.0;
    if (is_longitude) {
        double lat_rad = latitude * M_PI / 180.0;
        return degrees * METERS_PER_DEG_LAT * std::cos(lat_rad);
    } else {
        return degrees * METERS_PER_DEG_LAT;
    }
}
//这是一个工具函数，用于将经度或纬度（度）转换为米，或者将米转换为经度或纬度（度）。
double FixedHeightComparator::MetersToDegrees(double meters, double latitude, bool is_longitude) {
    const double METERS_PER_DEG_LAT = 111320.0;
    if (is_longitude) {
        double lat_rad = latitude * M_PI / 180.0;
        return meters / (METERS_PER_DEG_LAT * std::cos(lat_rad));
    } else {
        return meters / METERS_PER_DEG_LAT;
    }
}

// 【注意】CalculatePowerDensity函数已被删除
// 现在统一使用Radar类的CalculateSinglePointPowerDensityFromVec3d接口

FixedHeightComparator::ExperimentResult FixedHeightComparator::RunExperiment(
    Radar* radar,
    const Vec3d& radar_pos,
    double lon_range[2],
    double lat_range[2],
    double fixed_height,
    double grid_resolution_m,
    RTree3d* rtree,
    DEMLoader* dem
) {
    using namespace std::chrono;
    auto experiment_start = high_resolution_clock::now();

    std::cout << "\n========== Fixed-Height Plane Experiment (OpenMP Parallel) ==========" << std::endl;
    std::cout << "Fixed height: " << fixed_height << " m" << std::endl;
    std::cout << "Grid resolution: " << grid_resolution_m << " m" << std::endl;
    // 获取OpenMP线程数
    int num_threads = omp_get_max_threads();
    std::cout << "OpenMP threads: " << num_threads << "======================" <<std::endl;
    int thread = 4;
    omp_set_num_threads(thread);  // 设置为你想要的线程数，例如8
    ExperimentResult result;
    result.fixed_height = fixed_height;
    result.grid_resolution = grid_resolution_m;
    result.dem_file = ""; // 将在调用时设置
    result.dem_resolution = dem->GetResolution();
    result.num_threads = num_threads;

    // 计算网格参数
    double center_lat = (lat_range[0] + lat_range[1]) / 2.0;
    double lon_step = MetersToDegrees(grid_resolution_m, center_lat, true);
    double lat_step = MetersToDegrees(grid_resolution_m, center_lat, false);

    int steps_lon = static_cast<int>(std::ceil((lon_range[1] - lon_range[0]) / lon_step));
    int steps_lat = static_cast<int>(std::ceil((lat_range[1] - lat_range[0]) / lat_step));
    size_t total_points = static_cast<size_t>(steps_lon) * steps_lat;

    std::cout << "Grid dimensions: " << steps_lon << " x " << steps_lat
              << " = " << total_points << " points" << std::endl;

    // 预分配结果数组
    result.grid_points.resize(total_points);
    result.total_points = total_points;

    // RTree坐标系常量
    const double indexRange_x = -1975;
    const double indexRange_y = 43;
    const int min_x = 835000;
    const int min_y = 815500;

    // 主线程计算雷达局部坐标（只需计算一次）
    PJ_CONTEXT *main_proj_context = proj_context_create();
    PJ *main_proj = proj_create_crs_to_crs(main_proj_context, "EPSG:4326", "EPSG:2326", NULL);
    PJ_COORD radar_coord = proj_coord(radar_pos.y, radar_pos.x, 0, 0);
    PJ_COORD radar_result = proj_trans(main_proj, PJ_FWD, radar_coord);
    double radar_local_x = radar_result.xy.y - min_x;
    double radar_local_y = radar_result.xy.x - min_y;
    proj_destroy(main_proj);
    proj_context_destroy(main_proj_context);

    std::cout << "Running parallel computation (visibility + power density)..." << std::endl;

    // 统计变量（使用reduction）
    size_t local_visible_mesh = 0;
    size_t local_visible_dem = 0;
    size_t local_disagreement = 0;

    // min/max需要用critical section保护
    float global_min_density = 1e12f;
    float global_max_density = 0.0f;

    omp_set_num_threads(thread);  // 设置为你想要的线程数，例如8
    auto parallel_start = high_resolution_clock::now();  // 现在才开始计时

    // OpenMP并行循环：每个线程同时完成可见性判断和能量计算
    #pragma omp parallel reduction(+:local_visible_mesh, local_visible_dem, local_disagreement)
    {
    

        // 每个线程创建自己的其他资源（避免竞争）
        PJ_CONTEXT *thread_proj_context = proj_context_create();
        PJ *thread_proj = proj_create_crs_to_crs(thread_proj_context, "EPSG:4326", "EPSG:2326", NULL);
        caltools ct;
        DEMVisibility thread_dem_visibility(dem);

        // 线程局部的min/max
        float thread_min = 1e12f;
        float thread_max = 0.0f;

        #pragma omp for schedule(dynamic, 100)
        for (size_t idx = 0; idx < total_points; ++idx) {
            int i = idx / steps_lon;  // 行索引
            int j = idx % steps_lon;  // 列索引

            GridPoint& point = result.grid_points[idx];

            // 计算点的位置（EPSG:4326）
            point.position.x = lon_range[0] + j * lon_step;
            point.position.y = lat_range[0] + i * lat_step;
            point.position.z = fixed_height;

            // 转换目标点到EPSG:2326局部坐标
            PJ_COORD target_coord = proj_coord(point.position.y, point.position.x, 0, 0);
            PJ_COORD target_result = proj_trans(thread_proj, PJ_FWD, target_coord);
            double target_local_x = target_result.xy.y - min_x;
            double target_local_y = target_result.xy.x - min_y;

            // 1. 三角面片可见性判断（使用线程局部RTree，无需同步）
            Rect3d search_rect(
                target_local_x + indexRange_x,
                indexRange_y + target_local_y,
                fixed_height,
                indexRange_x + radar_local_x,
                indexRange_y + radar_local_y,
                radar->_cached_abs_alt_m
            );

            bool occluded_mesh = rtree->Intersect3d(search_rect.min, search_rect.max, IntersectCallback);
            point.visible_mesh = !occluded_mesh;

            // // 2. DEM可见性判断
            // bool occluded_dem = thread_dem_visibility.IsOccluded(radar_pos, point.position);
            // point.visible_dem = !occluded_dem;

            // 3. 能量密度计算（只对mesh可见的点计算）
            if (point.visible_mesh) {
                point.power_density = radar->CalculateSinglePointPowerDensity(point.position, &ct);

                // if (point.power_density > 0) {
                //     if (point.power_density < thread_min) thread_min = point.power_density;
                //     if (point.power_density > thread_max) thread_max = point.power_density;
                // }
            } else {
                point.power_density = 0.0f;
            }

            // // 统计
            // if (point.visible_mesh) local_visible_mesh++;
            // if (point.visible_dem) local_visible_dem++;
            // if (point.visible_mesh != point.visible_dem) local_disagreement++;
        }
        // 清理线程局部资源
        proj_destroy(thread_proj);
        proj_context_destroy(thread_proj_context);
        // 注意：thread_rtree不在这里删除，而是在parallel区域结束后统一清理

        // // 更新全局min/max（线程安全）
        // #pragma omp critical
        // {
        //     if (thread_min < global_min_density) global_min_density = thread_min;
        //     if (thread_max > global_max_density) global_max_density = thread_max;
        // }
    }
    auto parallel_end = high_resolution_clock::now();
    // 保存统计结果
    result.visible_mesh_count = local_visible_mesh;
    result.visible_dem_count = local_visible_dem;
    result.disagreement_count = local_disagreement;
    result.min_power_density = global_min_density;
    result.max_power_density = global_max_density;

    // 如果没有可见点,设置合理的范围
    if (result.min_power_density >= 1e12f) {
        result.min_power_density = 1e-6f;
        result.max_power_density = 1e-5f;
    }

    auto experiment_end = high_resolution_clock::now();

    // 计算时间（这里把所有阶段合并为一个）
    // result.total_time = duration_cast<duration<double>>(experiment_end - experiment_start).count();
    double parallel_time = duration_cast<duration<double>>(parallel_end - parallel_start).count();
    cout << "实验所用时长为" << parallel_time << "====================";
    // // 为了兼容性，将时间分配到各个阶段（实际是同时进行的）
    // result.mesh_visibility_time = parallel_time / 3.0;
    // result.dem_visibility_time = parallel_time / 3.0;
    // result.power_calculation_time = parallel_time / 3.0;

    // std::cout << "\n=== Experiment Completed ===" << std::endl;
    // std::cout << "Parallel computation time: " << parallel_time << " seconds" << std::endl;
    // std::cout << "Total time: " << result.total_time << " seconds" << std::endl;
    // std::cout << "Visible (Mesh): " << result.visible_mesh_count << " points" << std::endl;
    // std::cout << "Visible (DEM): " << result.visible_dem_count << " points" << std::endl;
    // std::cout << "Disagreements: " << result.disagreement_count << " points" << std::endl;

    return result;
}

void FixedHeightComparator::GenerateVisualization(
    const ExperimentResult& result,
    const std::string& output_file
) {
    std::cout << "Generating PNG visualization: " << output_file << std::endl;

    // 推断网格维度
    int steps_lon = 1;
    int steps_lat = 1;
    for (size_t i = 1; i < result.grid_points.size(); ++i) {
        if (result.grid_points[i].position.y == result.grid_points[0].position.y) {
            steps_lon++;
        } else {
            break;
        }
    }
    steps_lat = result.grid_points.size() / steps_lon;

    std::cout << "  Image dimensions: " << steps_lon << " x " << steps_lat << std::endl;

    // 创建PNG图像
    png_bytep* rowPointers = (png_bytep*)malloc(steps_lat * sizeof(png_bytep));
    for (int i = 0; i < steps_lat; i++) {
        rowPointers[i] = (png_bytep)malloc(steps_lon * 4 * sizeof(png_byte));
    }

    // 填充像素
    for (int i = 0; i < steps_lat; i++) {
        for (int j = 0; j < steps_lon; j++) {
            size_t index = i * steps_lon + j;
            const auto& point = result.grid_points[index];

            Color color;

            // 颜色编码逻辑
            if (point.visible_mesh && point.visible_dem) {
                // 两种方法都可见：显示能量密度
                color = mapDensityToColorLog(
                    point.power_density,
                    result.min_power_density,
                    result.max_power_density
                );
            } else if (!point.visible_mesh && !point.visible_dem) {
                // 两种方法都不可见：灰色
                color = {128, 128, 128};
            } else if (!point.visible_mesh && point.visible_dem) {
                // 仅DEM可见（假阴性）：蓝色
                color = {0, 0, 255};
            } else {
                // 仅mesh可见（假阳性）：绿色
                color = {0, 255, 0};
            }

            int ii = steps_lat - 1 - i; // 垂直翻转
            rowPointers[ii][j * 4 + 0] = color.r;
            rowPointers[ii][j * 4 + 1] = color.g;
            rowPointers[ii][j * 4 + 2] = color.b;
            rowPointers[ii][j * 4 + 3] = 255; // Alpha
        }
                    // 颜色编码逻辑：只基于DEM可见性
                //     if (point.visible_dem) {
                //         // DEM可见：显示能量密度
                //         color = mapDensityToColorLog(
                //             point.power_density,
                //             result.min_power_density,
                //             result.max_power_density
                //         );
                //     } else {
                //         // DEM不可见：灰色
                //         color = {128, 128, 128};
                //     }
        
                //     int ii = steps_lat - 1 - i; // 垂直翻转
                //     rowPointers[ii][j * 4 + 0] = color.r;
                //     rowPointers[ii][j * 4 + 1] = color.g;
                //     rowPointers[ii][j * 4 + 2] = color.b;
                //     rowPointers[ii][j * 4 + 3] = 255; // Alpha
                // }


    }

    // ============ 添加色标条（Colorbar） ============
    const int colorbar_width = 60;      // 色标条宽度
    const int colorbar_margin_lr = 15;  // 左右边距
    const int colorbar_margin_tb = 40;  // 上下边距
    const int total_width = steps_lon + colorbar_width + colorbar_margin_lr * 2;

    // 创建扩展的图像数组（包含色标条）
    png_bytep* extended_rows = (png_bytep*)malloc(steps_lat * sizeof(png_bytep));
    for (int i = 0; i < steps_lat; i++) {
        extended_rows[i] = (png_bytep)malloc(total_width * 4 * sizeof(png_byte));

        // 复制原始图像
        memcpy(extended_rows[i], rowPointers[i], steps_lon * 4 * sizeof(png_byte));

        // 填充右侧空白区域为白色
        for (int j = steps_lon; j < total_width; j++) {
            extended_rows[i][j * 4 + 0] = 255; // R
            extended_rows[i][j * 4 + 1] = 255; // G
            extended_rows[i][j * 4 + 2] = 255; // B
            extended_rows[i][j * 4 + 3] = 255; // A
        }
    }
    // 保存PNG
    FILE *fp = fopen(output_file.c_str(), "wb");
    if (!fp) {
        std::cerr << "Error: Could not open file " << output_file << std::endl;
        return;
    }

    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    png_infop info_ptr = png_create_info_struct(png_ptr);
    png_init_io(png_ptr, fp);
    png_set_IHDR(png_ptr, info_ptr, steps_lon, steps_lat, 8, PNG_COLOR_TYPE_RGB_ALPHA,
                 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    png_write_info(png_ptr, info_ptr);
    png_write_image(png_ptr, rowPointers);
    png_write_end(png_ptr, NULL);
    fclose(fp);
    png_destroy_write_struct(&png_ptr, &info_ptr);

    // 清理
    for (int i = 0; i < steps_lat; i++) {
        free(rowPointers[i]);
    }
    free(rowPointers);

    std::cout << "  PNG saved successfully" << std::endl;
}

void FixedHeightComparator::ExportResultsToCSV(
    const ExperimentResult& result,
    const std::string& output_file
) {
    std::cout << "Exporting results to CSV: " << output_file << std::endl;

    std::ofstream file(output_file);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << output_file << std::endl;
        return;
    }

    // 设置输出精度，确保坐标不被截断
    file << std::fixed << std::setprecision(10);

    // 写入CSV头
    file << "lon,lat,alt,visible_mesh,visible_dem,disagreement,power_density\n";

    // 写入数据
    for (const auto& point : result.grid_points) {
        file << point.position.x << ","
             << point.position.y << ","
             << point.position.z << ","
             << (point.visible_mesh ? 1 : 0) << ","
             << (point.visible_dem ? 1 : 0) << ","
             << (point.visible_mesh != point.visible_dem ? 1 : 0) << ","
             << point.power_density << "\n";
    }

    file.close();
    std::cout << "  CSV exported successfully" << std::endl;
}

void FixedHeightComparator::PrintSummary(const ExperimentResult& result) {
    std::cout << "\n========== Experiment Summary ==========" << std::endl;
    std::cout << "Fixed height: " << result.fixed_height << " m" << std::endl;
    std::cout << "Grid resolution: " << result.grid_resolution << " m" << std::endl;
    std::cout << "DEM file: " << result.dem_file << std::endl;
    std::cout << "DEM resolution: " << result.dem_resolution << " deg/pixel" << std::endl;
    std::cout << "\nVisibility Statistics:" << std::endl;
    std::cout << "  Total points: " << result.total_points << std::endl;
    std::cout << "  Visible (Mesh): " << result.visible_mesh_count
              << " (" << (100.0 * result.visible_mesh_count / result.total_points) << "%)" << std::endl;
    std::cout << "  Visible (DEM):  " << result.visible_dem_count
              << " (" << (100.0 * result.visible_dem_count / result.total_points) << "%)" << std::endl;
    std::cout << "  Disagreements:  " << result.disagreement_count
              << " (" << (100.0 * result.disagreement_count / result.total_points) << "%)" << std::endl;

    // 计算准确率
    double accuracy = 100.0 * (result.total_points - result.disagreement_count) / result.total_points;
    std::cout << "  DEM Accuracy:   " << accuracy << "%" << std::endl;

    std::cout << "\nPower Density Range:" << std::endl;
    std::cout << "  Min: " << result.min_power_density << " W/m²" << std::endl;
    std::cout << "  Max: " << result.max_power_density << " W/m²" << std::endl;
    std::cout << "========================================\n" << std::endl;
}

} // namespace hiradar
