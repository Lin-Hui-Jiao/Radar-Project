#include "hiradar/fixed_height_comparator.hpp"
#include "hiradar/dem_generator.hpp"
#include "hiradar/dem_loader.hpp"
#include "hiradar/RTree.h"
#include "hiradar/radar_pool.hpp"
#include "hiradar/grid.hpp"
#include "hiradar/occlusion_utils.h"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <proj.h>

using namespace hiradar;

namespace {

struct DEMConfig {
    double resolution_degree; //角度分辨率
    double approx_resolution_m; //近似分辨率（米）
    std::string description; //描述
};

const DEMConfig DEM_CONFIGS[] = {
    //{0.00001, 1.0, "Medium-precision DEM - UAV photogrammetry"},
    {0.00005, 5.0, "Low-precision DEM - SRTM/ASTER GDEM"},
    {0.0001, 10.0, "Coarse DEM - Digitized topographic map"},
    {0.00015, 15.0, "15m Coarse DEM - Digitized topographic map"},
    {0.0002, 20.0, "20m Coarse DEM - Digitized topographic map"},
};

struct ExperimentOptions {
    DatasetBounds bounds;
    RTreeCoordConfig coord_config;

    std::string index_file = "../../test_area_1_1.3idx";
    std::string output_dir = "fixed_height_2km_outputs";

    double radar_lon = 0.0;
    double radar_lat = 0.0;
    double radar_alt_relative_m = 80.0;
    double base_phi = 0.0;
    double base_theta = 0.0;
    double fixed_height_m = 50.0;

    bool radar_lon_set = false;
    bool radar_lat_set = false;
    bool show_help = false;

    ExperimentOptions() {
        radar_lon = (bounds.min_lon + bounds.max_lon) / 2.0;
        radar_lat = (bounds.min_lat + bounds.max_lat) / 2.0;
    } //  默认初始化雷达位置为区域中心
};

void PrintUsage(const char* exe_name) {
    std::cout
        << "Usage: " << exe_name << " [options]\n\n"
        << "Dataset options:\n"
        << "  --index <path>             RTree index file (default: ../../test_area_1_1.3idx)\n"
        << "  --output-dir <dir>         Output/cache directory (default: fixed_height_2km_outputs)\n"
        << "  --min-lon <deg>            Study area min longitude\n"
        << "  --min-lat <deg>            Study area min latitude\n"
        << "  --max-lon <deg>            Study area max longitude\n"
        << "  --max-lat <deg>            Study area max latitude\n"
        << "  --index-range-x <value>    RTree X offset\n"
        << "  --index-range-y <value>    RTree Y offset\n"
        << "  --min-x <value>            EPSG:2326 min X used by the index\n"
        << "  --min-y <value>            EPSG:2326 min Y used by the index\n\n"
        << "Radar/experiment options:\n"
        << "  --radar-lon <deg>          Radar longitude (default: dataset center)\n"
        << "  --radar-lat <deg>          Radar latitude (default: dataset center)\n"
        << "  --radar-alt <m>            Radar height above terrain (default: 80)\n"
        << "  --base-phi <deg>           Radar base pitch angle (default: 0)\n"
        << "  --base-theta <deg>         Radar base azimuth angle (default: 0)\n"
        << "  --fixed-height <m>         Fixed plane height ASL (default: 50)\n"
        << "  --help                     Show this help\n";
}
// 把命令行里传进来的字符串转成 double，并检查是否合法。
double ParseDouble(const char* value, const std::string& option_name) {
    char* end = nullptr;
    const double parsed = std::strtod(value, &end);
    if (end == value || *end != '\0') {
        throw std::invalid_argument("Invalid numeric value for " + option_name + ": " + value);
    }
    return parsed;
}

bool ParseOptions(int argc, char** argv, ExperimentOptions& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];

        auto require_value = [&](const std::string& option_name) -> const char* {
            if (i + 1 >= argc) {
                throw std::invalid_argument("Missing value for " + option_name);
            }
            return argv[++i];
        };

        if (arg == "--help" || arg == "-h") {
            options.show_help = true;
            return true;
        } else if (arg == "--index") {
            options.index_file = require_value(arg);
        } else if (arg == "--output-dir") {
            options.output_dir = require_value(arg);
        } else if (arg == "--min-lon") {
            options.bounds.min_lon = ParseDouble(require_value(arg), arg);
        } else if (arg == "--min-lat") {
            options.bounds.min_lat = ParseDouble(require_value(arg), arg);
        } else if (arg == "--max-lon") {
            options.bounds.max_lon = ParseDouble(require_value(arg), arg);
        } else if (arg == "--max-lat") {
            options.bounds.max_lat = ParseDouble(require_value(arg), arg);
        } else if (arg == "--index-range-x") {
            options.coord_config.index_range_x = ParseDouble(require_value(arg), arg);
        } else if (arg == "--index-range-y") {
            options.coord_config.index_range_y = ParseDouble(require_value(arg), arg);
        } else if (arg == "--min-x") {
            options.coord_config.min_x = ParseDouble(require_value(arg), arg);
        } else if (arg == "--min-y") {
            options.coord_config.min_y = ParseDouble(require_value(arg), arg);
        } else if (arg == "--radar-lon") {
            options.radar_lon = ParseDouble(require_value(arg), arg);
            options.radar_lon_set = true;
        } else if (arg == "--radar-lat") {
            options.radar_lat = ParseDouble(require_value(arg), arg);
            options.radar_lat_set = true;
        } else if (arg == "--radar-alt") {
            options.radar_alt_relative_m = ParseDouble(require_value(arg), arg);
        } else if (arg == "--base-phi") {
            options.base_phi = ParseDouble(require_value(arg), arg);
        } else if (arg == "--base-theta") {
            options.base_theta = ParseDouble(require_value(arg), arg);
        } else if (arg == "--fixed-height") {
            options.fixed_height_m = ParseDouble(require_value(arg), arg);
        } else {
            throw std::invalid_argument("Unknown option: " + arg);
        }
    }

    if (!options.radar_lon_set) {
        options.radar_lon = (options.bounds.min_lon + options.bounds.max_lon) / 2.0;
    }
    if (!options.radar_lat_set) {
        options.radar_lat = (options.bounds.min_lat + options.bounds.max_lat) / 2.0;
    }

    if (options.bounds.min_lon >= options.bounds.max_lon || options.bounds.min_lat >= options.bounds.max_lat) {
        throw std::invalid_argument("Invalid dataset bounds: min must be smaller than max");
    }

    return true;
}
// ResolutionTag(double resolution_m)
// 作用：把 DEM 分辨率转成文件名里用的字符串标签。
std::string ResolutionTag(double resolution_m) {
    const double rounded = std::round(resolution_m);
    if (std::abs(resolution_m - rounded) < 1e-6) {
        return std::to_string(static_cast<int>(rounded)) + "m";
    }

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(1) << resolution_m << "m";
    return oss.str();
}
//把经纬度坐标转换成 RTree 索引内部使用的坐标。4326->局部坐标系
bool LonLatToRTreeCoord(
    double lon,
    double lat,
    const RTreeCoordConfig& coord_config,
    double& rtree_x,
    double& rtree_y
) {
    PJ_CONTEXT* context = proj_context_create();
    PJ* proj_from = proj_create_crs_to_crs(context, "EPSG:4326", "EPSG:2326", nullptr); // 每次都要创建新的为什么不直接复用呢？作为参数传进来
    if (!proj_from) {
        proj_context_destroy(context);
        return false;
    }

    PJ_COORD coord = proj_coord(lat, lon, 0, 0);
    PJ_COORD result = proj_trans(proj_from, PJ_FWD, coord);

    proj_destroy(proj_from);
    proj_context_destroy(context);

    rtree_x = coord_config.index_range_x + (result.xy.y - coord_config.min_x);
    rtree_y = coord_config.index_range_y + (result.xy.x - coord_config.min_y);
    return true;
}
// 查询高度
bool QueryTerrainHeight(
    RTree3d* rtree,
    double lon,
    double lat,
    const RTreeCoordConfig& coord_config,
    double& terrain_height
) {
    double rtree_x = 0.0;
    double rtree_y = 0.0;
    if (!LonLatToRTreeCoord(lon, lat, coord_config, rtree_x, rtree_y)) {
        return false;
    }

    terrain_height = rtree->Getheight3d(rtree_x, rtree_y, HeightCallback);
    return terrain_height >= MINH;
}

} // namespace

int main(int argc, char** argv) {
    ExperimentOptions options;
    try {
        ParseOptions(argc, argv, options);
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n\n";
        PrintUsage(argv[0]);
        return 1;
    }

    if (options.show_help) {
        PrintUsage(argv[0]);
        return 0;
    }

    std::cout << "========================================" << std::endl;
    std::cout << "  Fixed-Height Visibility Comparison" << std::endl;
    std::cout << "  3D Triangle Mesh vs 2.5D DEM" << std::endl;
    std::cout << "========================================\n" << std::endl;

    std::cout << "[1/6] Loading triangle mesh data..." << std::endl;

    RTree3d* rtree = new RTree3d();
    if (!rtree->Load(options.index_file.c_str())) {
        std::cerr << "Error: Cannot load R-Tree index file: " << options.index_file << std::endl;
        std::cerr << "Use --index <path> to point to test_area_1_1.3idx if you run from a different directory." << std::endl;
        delete rtree;
        return 1;
    }

    std::cout << "  Index file: " << options.index_file << std::endl;
    std::cout << "  Triangle mesh loaded successfully\n" << std::endl;

    std::cout << "[2/6] Creating radar object..." << std::endl;

    Position radar_init_pos;
    DecimalToDMS(radar_init_pos.lon, options.radar_lon);
    DecimalToDMS(radar_init_pos.lat, options.radar_lat);
    radar_init_pos.alt = options.radar_alt_relative_m;

    const float base_phi = static_cast<float>(options.base_phi);
    const float base_theta = static_cast<float>(options.base_theta);
    Radar* radar = CityGuardRadar(radar_init_pos, base_phi, base_theta);
    radar->BindRtree(rtree);

    // Radar::Move still contains legacy dataset offsets in radar.hpp. Keep this
    // fixed-height path on the new config by setting the cached ASL height here.
    double terrain_alt = 0.0;
    if (QueryTerrainHeight(rtree, options.radar_lon, options.radar_lat, options.coord_config, terrain_alt)) {
        radar->_cached_abs_alt_m = terrain_alt + options.radar_alt_relative_m;
    } else {
        std::cerr << "Warning: failed to query terrain height for radar with the 2km coord config; "
                  << "using radar relative height as absolute fallback." << std::endl;
        radar->_cached_abs_alt_m = options.radar_alt_relative_m;
    }

    Vec3d radar_pos;
    radar_pos.x = options.radar_lon;
    radar_pos.y = options.radar_lat;
    radar_pos.z = radar->_cached_abs_alt_m;

    std::cout << "  Radar type: CityGuardRadar (PHASED)" << std::endl;
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "  Radar position: lon=" << radar_pos.x
              << ", lat=" << radar_pos.y
              << ", relative_alt=" << options.radar_alt_relative_m << " m"
              << ", absolute_alt=" << radar_pos.z << " m" << std::endl;
    std::cout << "  Base phi/theta: " << options.base_phi << " / " << options.base_theta << std::endl;
    std::cout << "  Radar initialized and bound to R-Tree\n" << std::endl;

    std::cout << "[3/6] Defining study area..." << std::endl;

    double lon_range[2] = {options.bounds.min_lon, options.bounds.max_lon};
    double lat_range[2] = {options.bounds.min_lat, options.bounds.max_lat};

    const double center_lat = (lat_range[0] + lat_range[1]) / 2.0;
    const double meters_per_deg_lon = 111320.0 * std::cos(center_lat * M_PI / 180.0);
    const double meters_per_deg_lat = 110704.0;
    const double area_width_m = (lon_range[1] - lon_range[0]) * meters_per_deg_lon;
    const double area_height_m = (lat_range[1] - lat_range[0]) * meters_per_deg_lat;

    std::cout << "  Coordinate System: EPSG:4326 (WGS84)" << std::endl;
    std::cout << "  Lon range: [" << lon_range[0] << ", " << lon_range[1] << "] deg" << std::endl;
    std::cout << "  Lat range: [" << lat_range[0] << ", " << lat_range[1] << "] deg" << std::endl;
    std::cout << "  Approximate area size: " << area_width_m << " m x "
              << area_height_m << " m" << std::endl;
    std::cout << "  Study area defined\n" << std::endl;

    std::cout << "[4/6] Checking/generating multi-resolution DEMs..." << std::endl;

    const std::filesystem::path output_dir(options.output_dir);
    std::filesystem::create_directories(output_dir);

    std::vector<std::string> dem_files;
    const int num_configs = sizeof(DEM_CONFIGS) / sizeof(DEM_CONFIGS[0]);

    for (int i = 0; i < num_configs; i++) {
        const DEMConfig& config = DEM_CONFIGS[i];
        const std::string filename = "2km_dem_" + ResolutionTag(config.approx_resolution_m) + ".tif";
        const std::filesystem::path dem_path = output_dir / filename;
        const std::string dem_path_str = dem_path.string();

        if (std::filesystem::exists(dem_path)) {
            std::cout << "  Found existing: " << dem_path_str << std::endl;
            dem_files.push_back(dem_path_str);
        } else {
            std::cout << "  Generating " << dem_path_str << "..." << std::endl;
            const bool success = DEMGenerator::GenerateDEMFromTriangles(
                rtree,
                lon_range,
                lat_range,
                config.approx_resolution_m,
                dem_path_str.c_str(),
                options.coord_config
            );
            if (success) {
                dem_files.push_back(dem_path_str);
            } else {
                std::cerr << "    Warning: Failed to generate " << dem_path_str << std::endl;
            }
        }
    }
    std::cout << "  DEM files ready (" << dem_files.size() << " files)\n" << std::endl;
    std::cout << "[5/6] Setting up experiment parameters..." << std::endl;
    const double fixed_height = options.fixed_height_m;
    const double grid_resolution_m = 1.0; // TODO 之后还要进行修改
    const double lon_step = grid_resolution_m / meters_per_deg_lon;
    const double lat_step = grid_resolution_m / meters_per_deg_lat;
    const int expected_points = static_cast<int>(
        std::ceil((lon_range[1] - lon_range[0]) / lon_step) *
        std::ceil((lat_range[1] - lat_range[0]) / lat_step)
    );

    std::cout << "  Fixed height plane: " << fixed_height << " m (ASL)" << std::endl;
    std::cout << "  Grid resolution: " << grid_resolution_m << " m" << std::endl;
    std::cout << "  Expected grid points: ~" << expected_points << std::endl;
    std::cout << "  Output/cache directory: " << output_dir.string() << std::endl;
    std::cout << "  Experiment parameters set\n" << std::endl;

    std::cout << "[6/6] Running fixed-height comparison experiments...\n" << std::endl;

    for (size_t i = 0; i < dem_files.size(); i++) {
        const std::string& dem_file = dem_files[i];
        const DEMConfig& config = DEM_CONFIGS[i];
        const std::string resolution_tag = ResolutionTag(config.approx_resolution_m);

        std::cout << "========================================" << std::endl;
        std::cout << "Experiment " << (i + 1) << "/" << dem_files.size()
                  << ": " << config.description << std::endl;
        std::cout << "DEM file: " << dem_file << std::endl;
        std::cout << "Resolution: " << config.resolution_degree << " deg/pixel (~"
                  << config.approx_resolution_m << " m)" << std::endl;
        std::cout << "========================================" << std::endl;

        DEMLoader dem;
        if (!dem.Load(dem_file)) {
            std::cerr << "Error: Cannot load DEM file, skipping..." << std::endl;
            continue;
        }

        auto result = FixedHeightComparator::RunExperiment(
            radar,
            radar_pos,
            lon_range,
            lat_range,
            fixed_height,
            grid_resolution_m,
            rtree,
            &dem,
            options.coord_config
        );

        result.dem_file = dem_file;
        FixedHeightComparator::PrintSummary(result);

        const std::filesystem::path csv_file =
            output_dir / ("2km_fixed_height_result_" + resolution_tag + ".csv");
        FixedHeightComparator::ExportResultsToCSV(result, csv_file.string());

        const std::filesystem::path png_file =
            output_dir / ("2km_fixed_height_viz_" + resolution_tag + ".png");
        FixedHeightComparator::GenerateVisualization(result, png_file.string());

        std::cout << std::endl;
    }

    std::cout << "Cleaning up..." << std::endl;
    delete radar;
    delete rtree;

    std::cout << "\n========================================" << std::endl;
    std::cout << "  All Experiments Completed!" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "\nGenerated files are under: " << output_dir.string() << std::endl;
    std::cout << "  DEM files: 2km_dem_*.tif" << std::endl;
    std::cout << "  CSV results: 2km_fixed_height_result_*.csv" << std::endl;
    std::cout << "  PNG visualizations: 2km_fixed_height_viz_*.png" << std::endl;

    return 0;
}
