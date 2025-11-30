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

    std::cout << "\n========== Fixed-Height Plane Experiment (OpenMP Parallel - Optimized) ==========" << std::endl;
    std::cout << "Fixed height: " << fixed_height << " m" << std::endl;
    std::cout << "Grid resolution: " << grid_resolution_m << " m" << std::endl;

    // 获取可用线程数，让用户可以通过环境变量 OMP_NUM_THREADS 控制
    int num_threads = omp_get_max_threads();
    std::cout << "OpenMP threads available: " << num_threads << std::endl;

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

    // ========== 阶段1: 预计算所有坐标转换 ==========
    auto precompute_start = high_resolution_clock::now();

    // 主线程计算雷达局部坐标
    PJ_CONTEXT *main_proj_context = proj_context_create();
    PJ *main_proj = proj_create_crs_to_crs(main_proj_context, "EPSG:4326", "EPSG:2326", NULL);
    PJ_COORD radar_coord = proj_coord(radar_pos.y, radar_pos.x, 0, 0);
    PJ_COORD radar_result = proj_trans(main_proj, PJ_FWD, radar_coord);
    double radar_local_x = radar_result.xy.y - min_x;
    double radar_local_y = radar_result.xy.x - min_y;

    // 预计算网格点的经纬度坐标
    std::vector<double> lon_coords(steps_lon);
    std::vector<double> lat_coords(steps_lat);
    for (int j = 0; j < steps_lon; ++j) {
        lon_coords[j] = lon_range[0] + j * lon_step;
    }
    for (int i = 0; i < steps_lat; ++i) {
        lat_coords[i] = lat_range[0] + i * lat_step;
    }

    // 预计算所有网格点的局部坐标（关键优化：将PROJ转换移出主循环）
    std::vector<double> target_local_x_arr(total_points);
    std::vector<double> target_local_y_arr(total_points);

    #pragma omp parallel
    {
        PJ_CONTEXT *thread_proj_context = proj_context_create();
        PJ *thread_proj = proj_create_crs_to_crs(thread_proj_context, "EPSG:4326", "EPSG:2326", NULL);

        #pragma omp for schedule(static)
        for (size_t idx = 0; idx < total_points; ++idx) {
            int i = idx / steps_lon;
            int j = idx % steps_lon;

            double lon = lon_coords[j];
            double lat = lat_coords[i];

            // 保存位置信息到结果数组
            result.grid_points[idx].position.x = lon;
            result.grid_points[idx].position.y = lat;
            result.grid_points[idx].position.z = fixed_height;

            // 坐标转换并缓存结果
            PJ_COORD target_coord = proj_coord(lat, lon, 0, 0);
            PJ_COORD target_result = proj_trans(thread_proj, PJ_FWD, target_coord);
            target_local_x_arr[idx] = target_result.xy.y - min_x;
            target_local_y_arr[idx] = target_result.xy.x - min_y;
        }

        proj_destroy(thread_proj);
        proj_context_destroy(thread_proj_context);
    }

    proj_destroy(main_proj);
    proj_context_destroy(main_proj_context);

    auto precompute_end = high_resolution_clock::now();
    double precompute_time = duration_cast<duration<double>>(precompute_end - precompute_start).count();
    std::cout << "Coordinate precomputation time: " << precompute_time << " seconds" << std::endl;

    // 缓存雷达参数
    const double radar_cached_alt = radar->_cached_abs_alt_m;

    // ========== 阶段2: 预加载 RTree 池（每个线程一个副本）==========
    auto rtree_load_start = high_resolution_clock::now();

    std::vector<RTree3d*> rtree_pool(num_threads);
    const char* rtree_file = "../../test_area.3idx";

    std::cout << "Loading " << num_threads << " RTree copies for parallel processing..." << std::endl;

    // 并行加载 RTree 副本
    #pragma omp parallel for
    for (int t = 0; t < num_threads; ++t) {
        rtree_pool[t] = new RTree3d();
        rtree_pool[t]->Load(rtree_file);
    }

    auto rtree_load_end = high_resolution_clock::now();
    double rtree_load_time = duration_cast<duration<double>>(rtree_load_end - rtree_load_start).count();
    std::cout << "RTree pool loading time: " << rtree_load_time << " seconds" << std::endl;

    std::cout << "Running parallel computation (visibility + power density)..." << std::endl;

    // 统计变量
    size_t local_visible_mesh = 0;
    size_t local_visible_dem = 0;
    size_t local_disagreement = 0;
    float global_min_density = 1e12f;
    float global_max_density = 0.0f;

    // ========== 阶段3: 主计算循环（每个线程使用私有RTree）==========
    auto parallel_start = high_resolution_clock::now();

    #pragma omp parallel reduction(+:local_visible_mesh, local_visible_dem, local_disagreement) \
                         reduction(min:global_min_density) reduction(max:global_max_density)
    {
        // 获取当前线程ID，使用对应的RTree副本
        int tid = omp_get_thread_num();
        RTree3d* my_rtree = rtree_pool[tid];

        // 线程私有的计算工具
        caltools ct;

        // 使用guided调度平衡负载
        #pragma omp for schedule(guided, 64)
        for (size_t idx = 0; idx < total_points; ++idx) {
            GridPoint& point = result.grid_points[idx];

            // 使用预计算的局部坐标
            double target_local_x = target_local_x_arr[idx];
            double target_local_y = target_local_y_arr[idx];

            // 1. 三角面片可见性判断（使用线程私有的RTree）
            Rect3d search_rect(
                target_local_x + indexRange_x,
                indexRange_y + target_local_y,
                fixed_height,
                indexRange_x + radar_local_x,
                indexRange_y + radar_local_y,
                radar_cached_alt
            );

            bool occluded_mesh = my_rtree->Intersect3d(search_rect.min, search_rect.max, IntersectCallback);
            point.visible_mesh = !occluded_mesh;

            // 2. 能量密度计算（只对mesh可见的点计算）
            if (point.visible_mesh) {
                point.power_density = radar->CalculateSinglePointPowerDensity(point.position, &ct);
                local_visible_mesh++;
            } else {
                point.power_density = 0.0f;
            }
        }
    }
    auto parallel_end = high_resolution_clock::now();

    // 清理 RTree 池
    for (int t = 0; t < num_threads; ++t) {
        delete rtree_pool[t];
    }
    rtree_pool.clear();
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

    // 计算时间
    double parallel_time = duration_cast<duration<double>>(parallel_end - parallel_start).count();
    double total_time = duration_cast<duration<double>>(experiment_end - experiment_start).count();
    result.total_time = total_time;

    // 输出计时信息
    std::cout << "\n=== Timing Results ===" << std::endl;
    std::cout << "Parallel computation time: " << parallel_time << " seconds" << std::endl;
    std::cout << "Total experiment time: " << total_time << " seconds" << std::endl;
    std::cout << "Points per second: " << static_cast<double>(total_points) / parallel_time << std::endl;
    std::cout << "Threads used: " << num_threads << std::endl;

    // 估算理想加速比（假设完美并行）
    double estimated_serial_time = parallel_time * num_threads;
    std::cout << "Estimated speedup: " << estimated_serial_time / parallel_time << "x (ideal: " << num_threads << "x)" << std::endl;

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
