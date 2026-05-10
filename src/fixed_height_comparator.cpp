#include "hiradar/fixed_height_comparator.hpp"
#include "hiradar/radar.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <omp.h>
#include <png.h>
#include <proj.h>

namespace hiradar {

namespace {

struct Color {
    unsigned char r;
    unsigned char g;
    unsigned char b;
};

struct ColorPoint {
    float value;
    int r;
    int g;
    int b;
};

const std::vector<ColorPoint> kColormap = {
    {0.00f,  68,   1,  84},
    {0.10f,  71,  15, 107},
    {0.20f,  72,  36, 117},
    {0.30f,  67,  55, 125},
    {0.40f,  56,  75, 126},
    {0.50f,  42,  95, 125},
    {0.60f,  32, 115, 117},
    {0.70f,  34, 135,  97},
    {0.80f,  68, 154,  68},
    {0.90f, 122, 170,  40},
    {1.00f, 253, 231,  37}
};

PJ* CreateProjectedTransform(PJ_CONTEXT* context, const char* source_crs, const char* target_crs) {
    PJ* raw = proj_create_crs_to_crs(context, source_crs, target_crs, nullptr);
    if (!raw) {
        return nullptr;
    }

    PJ* normalized = proj_normalize_for_visualization(context, raw);
    if (normalized) {
        proj_destroy(raw);
        return normalized;
    }

    return raw;
}

bool LonLatToProjected(PJ* transform, double lon, double lat, double& x2326, double& y2326) {
    PJ_COORD coord = proj_coord(lon, lat, 0, 0);
    PJ_COORD result = proj_trans(transform, PJ_FWD, coord);
    x2326 = result.xy.x;
    y2326 = result.xy.y;
    return std::isfinite(x2326) && std::isfinite(y2326);
}

bool ProjectedToLonLat(PJ* transform, double x2326, double y2326, double& lon, double& lat) {
    PJ_COORD coord = proj_coord(x2326, y2326, 0, 0);
    PJ_COORD result = proj_trans(transform, PJ_FWD, coord);
    lon = result.xy.x;
    lat = result.xy.y;
    return std::isfinite(lon) && std::isfinite(lat);
}

Color Lerp(const ColorPoint& c1, const ColorPoint& c2, float value) {
    const float t = (value - c1.value) / (c2.value - c1.value);
    return {
        static_cast<unsigned char>(c1.r + t * (c2.r - c1.r)),
        static_cast<unsigned char>(c1.g + t * (c2.g - c1.g)),
        static_cast<unsigned char>(c1.b + t * (c2.b - c1.b))
    };
}

Color MapDensityToColorLog(float density, float min_density, float max_density) {
    if (density <= 0.0f || max_density <= min_density || max_density <= 0.0f || min_density <= 0.0f) {
        return {128, 128, 128};
    }

    const float log_density = std::log(density);
    const float log_min = std::log(min_density);
    const float log_max = std::log(max_density);
    if (std::abs(log_max - log_min) < 1e-6f) {
        const auto& mid = kColormap[kColormap.size() / 2];
        return {
            static_cast<unsigned char>(mid.r),
            static_cast<unsigned char>(mid.g),
            static_cast<unsigned char>(mid.b)
        };
    }

    float normalized_v = (log_density - log_min) / (log_max - log_min);
    normalized_v = std::max(0.0f, std::min(1.0f, normalized_v));
    normalized_v = std::pow(normalized_v, 0.6f);

    if (normalized_v <= kColormap.front().value) {
        const auto& front = kColormap.front();
        return {
            static_cast<unsigned char>(front.r),
            static_cast<unsigned char>(front.g),
            static_cast<unsigned char>(front.b)
        };
    }
    if (normalized_v >= kColormap.back().value) {
        const auto& back = kColormap.back();
        return {
            static_cast<unsigned char>(back.r),
            static_cast<unsigned char>(back.g),
            static_cast<unsigned char>(back.b)
        };
    }

    for (size_t i = 0; i + 1 < kColormap.size(); ++i) {
        if (normalized_v >= kColormap[i].value && normalized_v <= kColormap[i + 1].value) {
            return Lerp(kColormap[i], kColormap[i + 1], normalized_v);
        }
    }

    const auto& back = kColormap.back();
    return {
        static_cast<unsigned char>(back.r),
        static_cast<unsigned char>(back.g),
        static_cast<unsigned char>(back.b)
    };
}

} // namespace

FixedHeightComparator::ExperimentResult FixedHeightComparator::RunExperiment(
    Radar* radar,
    const Vec3d& radar_pos,
    double lon_range[2],
    double lat_range[2],
    double fixed_height,
    double grid_resolution_m,
    RTree3d* rtree,
    DEMLoader* dem,
    const RTreeCoordConfig& coord_config
) {
    (void)lon_range;
    (void)lat_range;

    using namespace std::chrono;
    const auto experiment_start = high_resolution_clock::now();

    std::cout << "\n========== Fixed-Height Plane Experiment (EPSG:2326 Meter Grid) ==========" << std::endl;
    std::cout << "Fixed height: " << fixed_height << " m" << std::endl;
    std::cout << "Grid resolution: " << grid_resolution_m << " m" << std::endl;

    const int num_threads = omp_get_max_threads();
    std::cout << "OpenMP threads available: " << num_threads << std::endl;

    ExperimentResult result{};
    result.fixed_height = fixed_height;
    result.grid_resolution = grid_resolution_m;
    result.dem_file = "";
    result.dem_resolution = dem ? dem->GetResolution() : 0.0;
    result.num_threads = num_threads;

    if (grid_resolution_m <= 0.0 || coord_config.area_width_m <= 0.0 || coord_config.area_height_m <= 0.0) {
        std::cerr << "Error: invalid fixed-height grid configuration" << std::endl;
        return result;
    }

    const int steps_x = static_cast<int>(std::ceil(coord_config.area_width_m / grid_resolution_m));
    const int steps_y = static_cast<int>(std::ceil(coord_config.area_height_m / grid_resolution_m));
    const size_t total_points = static_cast<size_t>(steps_x) * static_cast<size_t>(steps_y);

    result.grid_width = steps_x;
    result.grid_height = steps_y;
    result.total_points = total_points;

    std::cout << "Grid dimensions: " << steps_x << " x " << steps_y
              << " = " << total_points << " points" << std::endl;
    std::cout << "Projected X range: [" << coord_config.min_x << ", "
              << (coord_config.min_x + coord_config.area_width_m) << "] m" << std::endl;
    std::cout << "Projected Y range: [" << coord_config.min_y << ", "
              << (coord_config.min_y + coord_config.area_height_m) << "] m" << std::endl;

    if (steps_x <= 0 || steps_y <= 0 || total_points == 0) {
        std::cerr << "Error: invalid grid dimensions" << std::endl;
        return result;
    }

    result.grid_points.resize(total_points);

    const auto precompute_start = high_resolution_clock::now();

    PJ_CONTEXT* main_context = proj_context_create();
    PJ* lonlat_to_projected = CreateProjectedTransform(main_context, "EPSG:4326", "EPSG:2326");
    if (!lonlat_to_projected) {
        std::cerr << "Error: failed to create PROJ transformer (EPSG:4326 -> EPSG:2326)" << std::endl;
        proj_context_destroy(main_context);
        return result;
    }

    double radar_x2326 = 0.0;
    double radar_y2326 = 0.0;
    if (!LonLatToProjected(lonlat_to_projected, radar_pos.x, radar_pos.y, radar_x2326, radar_y2326)) {
        std::cerr << "Error: failed to project radar position to EPSG:2326" << std::endl;
        proj_destroy(lonlat_to_projected);
        proj_context_destroy(main_context);
        return result;
    }

    proj_destroy(lonlat_to_projected);
    proj_context_destroy(main_context);

    const double radar_rtree_x = coord_config.index_range_x + (radar_x2326 - coord_config.min_x);
    const double radar_rtree_y = coord_config.index_range_y + (radar_y2326 - coord_config.min_y);
    const double radar_cached_alt = radar ? radar->_cached_abs_alt_m : radar_pos.z;

    bool projection_failed = false;

    #pragma omp parallel
    {
        PJ_CONTEXT* thread_context = proj_context_create();
        PJ* projected_to_lonlat = CreateProjectedTransform(thread_context, "EPSG:2326", "EPSG:4326");
        bool thread_projection_failed = (projected_to_lonlat == nullptr);

        #pragma omp for schedule(static)
        for (size_t idx = 0; idx < total_points; ++idx) {
            const int row = static_cast<int>(idx / steps_x);
            const int col = static_cast<int>(idx % steps_x);

            const double x2326 = coord_config.min_x + (static_cast<double>(col) + 0.5) * grid_resolution_m;
            const double y2326 = coord_config.min_y + coord_config.area_height_m
                               - (static_cast<double>(row) + 0.5) * grid_resolution_m;

            GridPoint& point = result.grid_points[idx];
            point.x2326 = x2326;
            point.y2326 = y2326;
            point.rtree_x = coord_config.index_range_x + (x2326 - coord_config.min_x);
            point.rtree_y = coord_config.index_range_y + (y2326 - coord_config.min_y);
            point.position.z = fixed_height;
            point.visible_mesh = false;
            point.visible_dem = false;
            point.power_density = 0.0f;

            if (!thread_projection_failed) {
                double lon = 0.0;
                double lat = 0.0;
                if (ProjectedToLonLat(projected_to_lonlat, x2326, y2326, lon, lat)) {
                    point.position.x = lon;
                    point.position.y = lat;
                } else {
                    thread_projection_failed = true;
                }
            }
        }

        if (projected_to_lonlat) {
            proj_destroy(projected_to_lonlat);
        }
        proj_context_destroy(thread_context);

        if (thread_projection_failed) {
            #pragma omp critical
            {
                projection_failed = true;
            }
        }
    }

    if (projection_failed) {
        std::cerr << "Error: failed to convert one or more fixed-height points to EPSG:4326" << std::endl;
        return result;
    }

    const auto precompute_end = high_resolution_clock::now();
    const double precompute_time = duration_cast<duration<double>>(precompute_end - precompute_start).count();
    std::cout << "Coordinate precomputation time: " << precompute_time << " seconds" << std::endl;
    std::cout << "Using loaded RTree index for mesh visibility." << std::endl;

    DEMVisibility dem_visibility(dem);

    Vec3d radar_pos_for_dem;
    radar_pos_for_dem.x = radar_x2326;
    radar_pos_for_dem.y = radar_y2326;
    radar_pos_for_dem.z = radar_cached_alt;

    size_t local_visible_mesh = 0;
    size_t local_visible_dem = 0;
    size_t local_disagreement = 0;
    float global_min_density = std::numeric_limits<float>::max();
    float global_max_density = 0.0f;

    std::cout << "Running parallel computation (mesh visibility + DEM visibility + power density)..." << std::endl;
    const auto parallel_start = high_resolution_clock::now();

    #pragma omp parallel reduction(+:local_visible_mesh, local_visible_dem, local_disagreement) \
                         reduction(min:global_min_density) reduction(max:global_max_density)
    {
        caltools ct;

        #pragma omp for schedule(guided, 64)
        for (size_t idx = 0; idx < total_points; ++idx) {
            GridPoint& point = result.grid_points[idx];

            double start_point[3] = {
                radar_rtree_x,
                radar_rtree_y,
                radar_cached_alt
            };
            double end_point[3] = {
                point.rtree_x,
                point.rtree_y,
                fixed_height
            };

            const bool occluded_mesh = rtree
                ? rtree->Intersect3d(start_point, end_point, IntersectCallback)
                : true;
            point.visible_mesh = !occluded_mesh;

            Vec3d target_pos_for_dem;
            target_pos_for_dem.x = point.x2326;
            target_pos_for_dem.y = point.y2326;
            target_pos_for_dem.z = fixed_height;
            const bool occluded_dem = dem
                ? dem_visibility.IsOccludedProjected(radar_pos_for_dem, target_pos_for_dem)
                : true;
            point.visible_dem = !occluded_dem;

            if (point.visible_mesh) {
                local_visible_mesh++;
            }
            if (point.visible_dem) {
                local_visible_dem++;
            }
            if (point.visible_mesh != point.visible_dem) {
                local_disagreement++;
            }

            if (point.visible_mesh && radar) {
                point.power_density = radar->CalculateSinglePointPowerDensity(point.position, &ct);
                if (point.power_density > 0.0f) {
                    global_min_density = std::min(global_min_density, point.power_density);
                    global_max_density = std::max(global_max_density, point.power_density);
                }
            } else {
                point.power_density = 0.0f;
            }
        }
    }

    const auto parallel_end = high_resolution_clock::now();

    result.visible_mesh_count = local_visible_mesh;
    result.visible_dem_count = local_visible_dem;
    result.disagreement_count = local_disagreement;
    result.min_power_density = global_min_density;
    result.max_power_density = global_max_density;

    if (result.min_power_density == std::numeric_limits<float>::max()) {
        result.min_power_density = 1e-6f;
        result.max_power_density = 1e-5f;
    }

    const auto experiment_end = high_resolution_clock::now();
    const double parallel_time = duration_cast<duration<double>>(parallel_end - parallel_start).count();
    const double total_time = duration_cast<duration<double>>(experiment_end - experiment_start).count();

    result.total_time = total_time;
    result.mesh_visibility_time = 0.0;
    result.dem_visibility_time = 0.0;
    result.power_calculation_time = 0.0;

    std::cout << "\n=== Timing Results ===" << std::endl;
    std::cout << "Parallel computation time: " << parallel_time << " seconds" << std::endl;
    std::cout << "Total experiment time: " << total_time << " seconds" << std::endl;
    if (parallel_time > 0.0) {
        std::cout << "Points per second: " << static_cast<double>(total_points) / parallel_time << std::endl;
    }
    std::cout << "Threads used: " << num_threads << std::endl;

    return result;
}

void FixedHeightComparator::GenerateVisualization(
    const ExperimentResult& result,
    const std::string& output_file
) {
    std::cout << "Generating PNG visualization: " << output_file << std::endl;

    if (result.grid_points.empty() || result.grid_width <= 0 || result.grid_height <= 0) {
        std::cerr << "Error: no grid points available for visualization" << std::endl;
        return;
    }

    const int steps_x = result.grid_width;
    const int steps_y = result.grid_height;

    std::cout << "  Image dimensions: " << steps_x << " x " << steps_y << std::endl;

    std::vector<std::vector<png_byte>> rows(steps_y, std::vector<png_byte>(steps_x * 4));
    std::vector<png_bytep> row_pointers(steps_y);

    for (int i = 0; i < steps_y; i++) {
        row_pointers[i] = rows[i].data();
    }

    for (int row = 0; row < steps_y; row++) {
        for (int col = 0; col < steps_x; col++) {
            const size_t index = static_cast<size_t>(row) * steps_x + col;
            const auto& point = result.grid_points[index];

            Color color{};
            if (point.visible_mesh && point.visible_dem) {
                color = MapDensityToColorLog(
                    point.power_density,
                    result.min_power_density,
                    result.max_power_density
                );
            } else if (!point.visible_mesh && !point.visible_dem) {
                color = {128, 128, 128};
            } else if (!point.visible_mesh && point.visible_dem) {
                color = {0, 0, 255};
            } else {
                color = {0, 255, 0};
            }

            rows[row][col * 4 + 0] = color.r;
            rows[row][col * 4 + 1] = color.g;
            rows[row][col * 4 + 2] = color.b;
            rows[row][col * 4 + 3] = 255;
        }
    }

    FILE* fp = fopen(output_file.c_str(), "wb");
    if (!fp) {
        std::cerr << "Error: Could not open file " << output_file << std::endl;
        return;
    }

    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    if (!png_ptr) {
        fclose(fp);
        std::cerr << "Error: Could not create PNG write struct" << std::endl;
        return;
    }

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
        png_destroy_write_struct(&png_ptr, nullptr);
        fclose(fp);
        std::cerr << "Error: Could not create PNG info struct" << std::endl;
        return;
    }

    png_init_io(png_ptr, fp);
    png_set_IHDR(
        png_ptr, info_ptr, steps_x, steps_y, 8, PNG_COLOR_TYPE_RGB_ALPHA,
        PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE
    );
    png_write_info(png_ptr, info_ptr);
    png_write_image(png_ptr, row_pointers.data());
    png_write_end(png_ptr, nullptr);

    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);

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

    file << std::fixed << std::setprecision(10);
    file << "lon,lat,x2326,y2326,rtree_x,rtree_y,alt,visible_mesh,visible_dem,disagreement,power_density\n";

    for (const auto& point : result.grid_points) {
        file << point.position.x << ","
             << point.position.y << ","
             << point.x2326 << ","
             << point.y2326 << ","
             << point.rtree_x << ","
             << point.rtree_y << ","
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
    std::cout << "Grid dimensions: " << result.grid_width << " x " << result.grid_height << std::endl;
    std::cout << "DEM file: " << result.dem_file << std::endl;
    std::cout << "DEM resolution: " << result.dem_resolution << " m/pixel" << std::endl;
    std::cout << "\nVisibility Statistics:" << std::endl;
    std::cout << "  Total points: " << result.total_points << std::endl;

    const double denom = result.total_points > 0 ? static_cast<double>(result.total_points) : 1.0;
    std::cout << "  Visible (Mesh): " << result.visible_mesh_count
              << " (" << (100.0 * result.visible_mesh_count / denom) << "%)" << std::endl;
    std::cout << "  Visible (DEM):  " << result.visible_dem_count
              << " (" << (100.0 * result.visible_dem_count / denom) << "%)" << std::endl;
    std::cout << "  Disagreements:  " << result.disagreement_count
              << " (" << (100.0 * result.disagreement_count / denom) << "%)" << std::endl;

    const double accuracy = 100.0 * (denom - result.disagreement_count) / denom;
    std::cout << "  DEM Accuracy:   " << accuracy << "%" << std::endl;

    std::cout << "\nPower Density Range:" << std::endl;
    std::cout << "  Min: " << result.min_power_density << " W/m^2" << std::endl;
    std::cout << "  Max: " << result.max_power_density << " W/m^2" << std::endl;
    std::cout << "========================================\n" << std::endl;
}

} // namespace hiradar
