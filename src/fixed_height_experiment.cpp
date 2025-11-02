#include "hiradar/fixed_height_comparator.hpp"
#include "hiradar/dem_generator.hpp"
#include "hiradar/dem_loader.hpp"
#include "hiradar/RTree.h"
#include "hiradar/radar_pool.hpp"  // 用于创建Radar对象
#include "hiradar/grid.hpp"        // 用于DecimalToDMS函数
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

using namespace hiradar;

// DEM分辨率配置（与comparison_experiment.cpp保持一致）
struct DEMConfig {
    double resolution_degree;
    double approx_resolution_m;
    std::string description;
};

const DEMConfig DEM_CONFIGS[] = {
    // {0.000005, 0.5,  "High-precision DEM - LiDAR simulation"},
    {0.00001,  1.0,  "Medium-precision DEM - UAV photogrammetry"},
    // {0.00005,  5.0,  "Low-precision DEM - SRTM/ASTER GDEM"},
    // {0.0001,   10.0, "Coarse DEM - Digitized topographic map"},
    // {0.00015,  15.0, "15m Coarse DEM - Digitized topographic map"},
    // {0.0002,  20.0, "20m Coarse DEM - Digitized topographic map"},
};

int main(int argc, char** argv) {
    std::cout << "========================================" << std::endl;
    std::cout << "  Fixed-Height Visibility Comparison" << std::endl;
    std::cout << "  3D Triangle Mesh vs 2.5D DEM" << std::endl;
    std::cout << "========================================\n" << std::endl;

    // ========== 第1步：加载三角面片数据 ==========
    std::cout << "[1/5] Loading triangle mesh data..." << std::endl;

    RTree3d* rtree = new RTree3d();
    const char* rtree_file = "../../test_area.3idx";

    if (!rtree->Load(rtree_file)) {
        std::cerr << "Error: Cannot load R-Tree index file: " << rtree_file << std::endl;
        std::cerr << "Please ensure test_area.3idx exists in the current directory" << std::endl;
        delete rtree;
        return 1;
    }

    std::cout << "✓ Triangle mesh loaded successfully\n" << std::endl;

    // ========== 第2步：创建和初始化Radar对象 ==========
    std::cout << "[2/6] Creating radar object..." << std::endl;

    // 创建雷达位置（先创建临时位置用于初始化）
    Position radar_init_pos;
    DecimalToDMS(radar_init_pos.lon, 114.1670);  // 区域中心经度
    DecimalToDMS(radar_init_pos.lat, 22.2806);   // 区域中心纬度
    radar_init_pos.alt = 80.0;  // 雷达相对地面高度80米

    // 使用CityGuardRadar创建雷达对象
    float base_phi = 0.0f;   // 载具基础俯仰角
    float base_theta = 0.0f; // 载具基础方位角
    Radar* radar = CityGuardRadar(radar_init_pos, base_phi, base_theta);

    // 绑定R-Tree（用于遮挡检测）
    radar->BindRtree(rtree);

    std::cout << "  Radar type: CityGuardRadar (PHASED)" << std::endl;
    std::cout << "  Radar position: lon=" << 114.1670 << ", lat=" << 22.2806 << ", alt=80m" << std::endl;
    std::cout << "✓ Radar initialized and bound to R-Tree\n" << std::endl;

    // ========== 第3步：定义研究区域范围 ==========
    std::cout << "[3/6] Defining study area..." << std::endl;

    // 使用经纬度范围（EPSG:4326坐标系）
    double lon_range[2] = {114.164571, 114.169422};  // [min_lon, max_lon]
    double lat_range[2] = {22.278365, 22.28288};    // [min_lat, max_lat]

    std::cout << "  Coordinate System: EPSG:4326 (WGS84)" << std::endl;
    std::cout << "  Lon range: [" << lon_range[0] << ", " << lon_range[1] << "] deg" << std::endl;
    std::cout << "  Lat range: [" << lat_range[0] << ", " << lat_range[1] << "] deg" << std::endl;

    // 计算研究区域的近似尺寸
    double center_lat = (lat_range[0] + lat_range[1]) / 2.0;
    double meters_per_deg_lon = 111320.0 * std::cos(center_lat * M_PI / 180.0);
    double meters_per_deg_lat = 111320.0;
    // double area_width_m = (lon_range[1] - lon_range[0]) * meters_per_deg_lon;
    // double area_height_m = (lat_range[1] - lat_range[0]) * meters_per_deg_lat;
    //之后可能进行修改
    double area_width_m = 500;
    double area_height_m = 500;

    std::cout << "  Approximate area size: " << area_width_m << " m × "
              << area_height_m << " m" << std::endl;
    std::cout << "✓ Study area defined\n" << std::endl;

    // ========== 第4步：生成多分辨率DEM（如果尚未生成）==========
    std::cout << "[4/6] Checking/generating multi-resolution DEMs..." << std::endl;

    std::vector<std::string> dem_files;
    int num_configs = sizeof(DEM_CONFIGS) / sizeof(DEM_CONFIGS[0]);

    for (int i = 0; i < num_configs; i++) {
        const DEMConfig& config = DEM_CONFIGS[i];
        std::string filename = "dem_" +
                              std::to_string(static_cast<int>(config.approx_resolution_m)) + "m.tif";

        // 检查文件是否已存在
        FILE* test_file = fopen(filename.c_str(), "r");
        if (test_file) {
            fclose(test_file);
            std::cout << "  ✓ Found existing: " << filename << std::endl;
            dem_files.push_back(filename);
        } else {
            std::cout << "  Generating " << filename << "..." << std::endl;
            bool success = DEMGenerator::GenerateDEMFromTriangles(
                rtree, lon_range, lat_range, config.approx_resolution_m, filename.c_str()
            );
            if (success) {
                dem_files.push_back(filename);
            } else {
                std::cerr << "    Warning: Failed to generate " << filename << std::endl;
            }
        }
    }
    std::cout << "✓ DEM files ready (" << dem_files.size() << " files)\n" << std::endl;

    // ========== 第5步：设置实验参数 ==========
    std::cout << "[5/6] Setting up experiment parameters..." << std::endl;

    // 雷达位置（EPSG:4326坐标系：经度、纬度、高程）
    Vec3d radar_pos;
    // radar_pos.x = (lon_range[0] + lon_range[1]) / 2.0;  // 区域中心经度
    // radar_pos.y = (lat_range[0] + lat_range[1]) / 2.0;  // 区域中心纬度
    radar_pos.x = 114.1670;
    radar_pos.y = 22.2806;
    radar_pos.z = 164.2;  // 雷达高度（米，绝对海拔高度）

    // 重要：调用radar->Move()更新雷达位置（初始化缓存的坐标转换结果）
    Position radar_pos_dms;
    DecimalToDMS(radar_pos_dms.lon, radar_pos.x);
    DecimalToDMS(radar_pos_dms.lat, radar_pos.y);
    radar_pos_dms.alt = 80;
    radar->Move(radar_pos_dms, base_phi, base_theta, 0);

    std::cout << "  Radar position (EPSG:4326):" << std::endl;
    std::cout << "    Longitude: " << radar_pos.x << " deg" << std::endl;
    std::cout << "    Latitude:  " << radar_pos.y << " deg" << std::endl;
    std::cout << "    Altitude:  " << radar_pos.z << " m (absolute sea level)" << std::endl;

    // 固定高度平面参数
    double fixed_height = 50.0;  // 固定高度50米（海拔）
    std::cout << "请输入png图像的固定高度: ";
    std::cin >> fixed_height;
    double grid_resolution_m = 1.0;  // 网格分辨率1米

    std::cout << "  Fixed height plane: " << fixed_height << " m (above sea level)" << std::endl;
    std::cout << "  Grid resolution: " << grid_resolution_m << " m" << std::endl;

    // 计算预期网格大小
    double lon_step = grid_resolution_m / meters_per_deg_lon;
    double lat_step = grid_resolution_m / meters_per_deg_lat;
    int expected_points = static_cast<int>(
        std::ceil((lon_range[1] - lon_range[0]) / lon_step) *
        std::ceil((lat_range[1] - lat_range[0]) / lat_step)
    );
    std::cout << "  Expected grid points: ~" << expected_points << std::endl;

    std::cout << "✓ Experiment parameters set\n" << std::endl;

    // ========== 第6步：对每个DEM执行对比实验 ==========
    std::cout << "[6/6] Running fixed-height comparison experiments...\n" << std::endl;

    for (size_t i = 0; i < dem_files.size(); i++) {
        const std::string& dem_file = dem_files[i];
        const DEMConfig& config = DEM_CONFIGS[i];

        std::cout << "========================================" << std::endl;
        std::cout << "Experiment " << (i+1) << "/" << dem_files.size()
                  << ": " << config.description << std::endl;
        std::cout << "DEM file: " << dem_file << std::endl;
        std::cout << "Resolution: " << config.resolution_degree << " deg/pixel (~"
                  << config.approx_resolution_m << " m)" << std::endl;
        std::cout << "========================================" << std::endl;

        // 加载DEM
        DEMLoader dem;
        if (!dem.Load(dem_file)) {
            std::cerr << "Error: Cannot load DEM file, skipping..." << std::endl;
            continue;
        }

        // 执行固定高度对比实验（传入radar指针用于能量密度计算）
        auto result = FixedHeightComparator::RunExperiment(
            radar,          // 【新增】雷达对象指针
            radar_pos,
            lon_range,
            lat_range,
            fixed_height,
            grid_resolution_m,
            rtree,
            &dem
        );

        // 设置DEM文件名
        result.dem_file = dem_file;

        // 打印结果摘要
        FixedHeightComparator::PrintSummary(result);

        // 导出CSV结果
        std::string csv_file = "fixed_height_result_" +
                              std::to_string(static_cast<int>(config.approx_resolution_m)) + "m.csv";
        FixedHeightComparator::ExportResultsToCSV(result, csv_file);

        // 生成可视化PNG
        std::string png_file = "fixed_height_viz_" +
                              std::to_string(static_cast<int>(config.approx_resolution_m)) + "m.png";
        FixedHeightComparator::GenerateVisualization(result, png_file);

        std::cout << std::endl;
    }

    // ========== 清理资源 ==========
    std::cout << "Cleaning up..." << std::endl;
    delete radar;
    delete rtree;

    std::cout << "\n========================================" << std::endl;
    std::cout << "  All Experiments Completed!" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "\nGenerated files:" << std::endl;
    std::cout << "  DEM files: dem_*.tif" << std::endl;
    std::cout << "  CSV results: fixed_height_result_*.csv" << std::endl;
    std::cout << "  PNG visualizations: fixed_height_viz_*.png" << std::endl;
    std::cout << "\nColor encoding in PNG files:" << std::endl;
    std::cout << "  Gray:   Both methods agree (occluded)" << std::endl;
    std::cout << "  Yellow-Red gradient: Energy density (both visible)" << std::endl;
    std::cout << "  Blue:   Only DEM visible (false negative)" << std::endl;
    std::cout << "  Green:  Only Mesh visible (false positive)" << std::endl;
    std::cout << "\nNext steps:" << std::endl;
    std::cout << "  1. Compare PNG files across different DEM resolutions" << std::endl;
    std::cout << "  2. Analyze CSV files for quantitative statistics" << std::endl;
    std::cout << "  3. Focus on green/blue regions (disagreements)" << std::endl;

    return 0;
}
