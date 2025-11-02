#include "hiradar/dem_generator.hpp"
#include "hiradar/dem_loader.hpp"
#include "hiradar/dem_visibility.hpp"
#include "hiradar/visibility_comparator.hpp"
#include "hiradar/RTree.h"
#include "hiradar/occlusion_utils.h"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

using namespace hiradar;

// DEM分辨率配置
// 注意：香港地区（约22°N）1度经纬度约等于111km
// 因此 0.00001度 ≈ 1.11米，0.000005度 ≈ 0.55米
struct DEMConfig {
    double resolution_degree;  // 单位：度/像素（EPSG:4326）
    double approx_resolution_m;  // 近似米数（用于文件命名）
    std::string description;
};

const DEMConfig DEM_CONFIGS[] = {
    {0.000005, 0.5,  "High-precision DEM - LiDAR simulation"},
    {0.00001,  1.0,  "Medium-precision DEM - UAV photogrammetry"},
    {0.00005,  5.0,  "Low-precision DEM - SRTM/ASTER GDEM"},
    {0.0001,   10.0, "Coarse DEM - Digitized topographic map"}
};

int main(int argc, char** argv) {
    std::cout << "========================================" << std::endl;
    std::cout << "  Visibility Algorithm Comparison" << std::endl;
    std::cout << "  3D Triangle Mesh vs 2.5D DEM" << std::endl;
    std::cout << "========================================\n" << std::endl;

    // ========== 第1步：加载三角面片数据 ==========
    std::cout << "[1/5] Loading triangle mesh data..." << std::endl;

    
    RTree3d* rtree = new RTree3d();
    const char* rtree_file = "test_area.3idx";

    if (!rtree->Load(rtree_file)) {
        std::cerr << "Error: Cannot load R-Tree index file: " << rtree_file << std::endl;
        std::cerr << "Please ensure test_area.3idx exists in the current directory" << std::endl;
        delete rtree;
        return 1;
    }

    std::cout << "✓ Triangle mesh loaded successfully\n" << std::endl;

    // ========== 第2步：定义研究区域范围 ==========
    std::cout << "[2/5] Defining study area..." << std::endl;

    // 使用经纬度范围（EPSG:4326坐标系）
    // 整个实验统一使用EPSG:4326坐标系，确保坐标系一致性
    double lon_range[2] = {114.15, 114.20};  // [min_lon, max_lon]（十进制度）
    double lat_range[2] = {22.25, 22.30};    // [min_lat, max_lat]（十进制度）

    std::cout << "  Coordinate System: EPSG:4326 (WGS84)" << std::endl;
    std::cout << "  Lon range: [" << lon_range[0] << ", " << lon_range[1] << "] deg" << std::endl;
    std::cout << "  Lat range: [" << lat_range[0] << ", " << lat_range[1] << "] deg" << std::endl;

    // 计算研究区域的近似尺寸（用于信息显示）
    double center_lat = (lat_range[0] + lat_range[1]) / 2.0;
    double meters_per_deg_lon = 111320.0 * std::cos(center_lat * M_PI / 180.0);
    double meters_per_deg_lat = 111320.0;
    double area_width_m = (lon_range[1] - lon_range[0]) * meters_per_deg_lon;
    double area_height_m = (lat_range[1] - lat_range[0]) * meters_per_deg_lat;

    std::cout << "  Approximate area size: " << area_width_m << " m × "
              << area_height_m << " m" << std::endl;
    std::cout << "✓ Study area defined\n" << std::endl;

    // ========== 第3步：生成多分辨率DEM ==========
    std::cout << "[3/5] Generating multi-resolution DEMs..." << std::endl;

    std::vector<std::string> dem_files;
    int num_configs = sizeof(DEM_CONFIGS) / sizeof(DEM_CONFIGS[0]);

    for (int i = 0; i < num_configs; i++) {
        const DEMConfig& config = DEM_CONFIGS[i];
        std::string filename = "dem_" +
                              std::to_string(static_cast<int>(config.approx_resolution_m)) + "m.tif";

        std::cout << "\n  [" << (i+1) << "/" << num_configs << "] Generating " << filename
                  << " (" << config.description << ")" << std::endl;
        std::cout << "    Resolution: " << config.resolution_degree << " deg/pixel (~"
                  << config.approx_resolution_m << " m)" << std::endl;

        // 使用新的EPSG:4326 DEM生成API
        bool success = DEMGenerator::GenerateDEMFromTriangles(
            rtree, lon_range, lat_range, config.resolution_degree, filename.c_str()
        );

        if (!success) {
            std::cerr << "Warning: Failed to generate " << filename << std::endl;
            continue;
        }

        dem_files.push_back(filename);
    }

    std::cout << "\n✓ DEM generation completed\n" << std::endl;

    // ========== 第4步：设置雷达位置和测试点 ==========
    std::cout << "[4/5] Setting up experiment parameters..." << std::endl;

    // 雷达位置（EPSG:4326坐标系：经度、纬度、高程）
    Vec3d radar_pos;
    radar_pos.x = (lon_range[0] + lon_range[1]) / 2.0;  // 区域中心经度（度）
    radar_pos.y = (lat_range[0] + lat_range[1]) / 2.0;  // 区域中心纬度（度）
    radar_pos.z = 80.0;  // 雷达高度（米，海拔高度）

    std::cout << "  Radar position (EPSG:4326):" << std::endl;
    std::cout << "    Longitude: " << radar_pos.x << " deg" << std::endl;
    std::cout << "    Latitude:  " << radar_pos.y << " deg" << std::endl;
    std::cout << "    Altitude:  " << radar_pos.z << " m" << std::endl;

    // 生成测试点
    int num_test_points = 10000;
    double max_range_m = 300.0;  // 300米探测范围
    double min_height = 1.0;     // 最小高度1米
    double max_height = 100.0;   // 最大高度100米

    // 将米距离转换为度距离（用于测试点生成）
    // 在给定纬度上，1度经度 ≈ 111320 * cos(lat) 米
    // 1度纬度 ≈ 111320 米
    double max_range_deg = max_range_m / meters_per_deg_lat;  // 使用纬度距离作为近似

    std::cout << "  Generating " << num_test_points << " test points..." << std::endl;
    std::cout << "    Max range: " << max_range_m << " m (~"
              << max_range_deg << " deg)" << std::endl;
    std::cout << "    Height range: [" << min_height << ", " << max_height << "] m" << std::endl;

    auto test_points = VisibilityComparator::GenerateTestPoints(
        radar_pos, max_range_deg, num_test_points, min_height, max_height, 42  // 固定种子42
    );

    std::cout << "✓ Experiment parameters set (" << test_points.size()
              << " points generated)\n" << std::endl;

    // ========== 第5步：对每个DEM执行对比实验 ==========
    std::cout << "[5/5] Running comparison experiments...\n" << std::endl;

    for (size_t i = 0; i < dem_files.size(); i++) {
        const std::string& dem_file = dem_files[i];
        const DEMConfig& config = DEM_CONFIGS[i];

        std::cout << "----------------------------------------" << std::endl;
        std::cout << "Experiment " << (i+1) << "/" << dem_files.size()
                  << ": " << config.description << std::endl;
        std::cout << "DEM file: " << dem_file << std::endl;
        std::cout << "Resolution: " << config.resolution_degree << " deg/pixel (~"
                  << config.approx_resolution_m << " m)" << std::endl;
        std::cout << "----------------------------------------" << std::endl;

        // 加载DEM
        DEMLoader dem;
        if (!dem.Load(dem_file)) {
            std::cerr << "Error: Cannot load DEM file, skipping..." << std::endl;
            continue;
        }

        // 执行对比实验
        auto result = VisibilityComparator::Compare(
            radar_pos, test_points, rtree, &dem
        );

        // 打印结果摘要
        VisibilityComparator::PrintSummary(result);

        // 导出CSV结果
        std::string csv_file = "result_" +
                              std::to_string(static_cast<int>(config.approx_resolution_m)) + "m.csv";
        VisibilityComparator::ExportResultsToCSV(result, csv_file);

        // 生成可视化（可选）
        std::string png_file = "disagreement_" +
                              std::to_string(static_cast<int>(config.approx_resolution_m)) + "m.png";
        VisibilityComparator::VisualizeDisagreements(result, radar_pos, png_file);

        std::cout << std::endl;
    }

    // ========== 清理资源 ==========
    std::cout << "\nCleaning up..." << std::endl;
    delete rtree;

    std::cout << "\n========================================" << std::endl;
    std::cout << "  Experiment Completed!" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "\nGenerated files:" << std::endl;
    std::cout << "  DEM files: dem_*.tif" << std::endl;
    std::cout << "  Result files: result_*.csv" << std::endl;
    std::cout << "  Visualization: disagreement_*.png (if implemented)" << std::endl;
    std::cout << "\nNext steps:" << std::endl;
    std::cout << "  1. Analyze CSV files for detailed statistics" << std::endl;
    std::cout << "  2. Create graphs and tables for your paper" << std::endl;
    std::cout << "  3. Visualize disagreement points on maps" << std::endl;

    return 0;
}
