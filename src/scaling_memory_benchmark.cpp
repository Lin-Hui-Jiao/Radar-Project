#include "hiradar/dem_loader.hpp"
#include "hiradar/grid.hpp"
#include "hiradar/occlusion_utils.h"
#include "hiradar/radar_pool.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <tuple>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <proj.h>

namespace {

struct Bounds {
    double min_lon = 114.164571;
    double max_lon = 114.169422;
    double min_lat = 22.278365;
    double max_lat = 22.282880;
};

struct RTreeCoordConfig {
    double index_range_x = -1975.0;
    double index_range_y = 43.0;
    double min_x = 835000.0;
    double min_y = 815500.0;
};

struct Options {
    Bounds bounds;
    RTreeCoordConfig coord;

    std::string method = "both";
    std::string index_file = "../../test_area.3idx";
    std::string dem_file = "dem_1m.tif";
    std::string csv_file = "scaling_memory_benchmark.csv";
    std::vector<std::size_t> sizes = {250000, 1000000, 5000000, 10000000};

    int threads = 0;
    int warmup_runs = 1;

    double min_alt_m = 0.0;
    double max_alt_m = 100.0;
    std::uint64_t seed = 20260517;
    double radar_lon = 114.1670;
    double radar_lat = 22.2806;
    double radar_alt_relative_m = 80.0;

    bool write_csv = true;
    bool show_help = false;
};

struct Projector {
    PJ_CONTEXT* context = nullptr;
    PJ* lonlat_to_projected = nullptr;

    Projector() {
        context = proj_context_create();
        PJ* raw = proj_create_crs_to_crs(context, "EPSG:4326", "EPSG:2326", nullptr);
        if (!raw) {
            return;
        }
        lonlat_to_projected = proj_normalize_for_visualization(context, raw);
        if (lonlat_to_projected) {
            proj_destroy(raw);
        } else {
            lonlat_to_projected = raw;
        }
    }

    ~Projector() {
        if (lonlat_to_projected) {
            proj_destroy(lonlat_to_projected);
        }
        if (context) {
            proj_context_destroy(context);
        }
    }

    bool IsValid() const {
        return context != nullptr && lonlat_to_projected != nullptr;
    }

    bool LonLatToProjected(double lon, double lat, double& x2326, double& y2326) const {
        if (!IsValid()) {
            return false;
        }
        PJ_COORD coord = proj_coord(lon, lat, 0, 0);
        PJ_COORD result = proj_trans(lonlat_to_projected, PJ_FWD, coord);
        x2326 = result.xy.x;
        y2326 = result.xy.y;
        return std::isfinite(x2326) && std::isfinite(y2326);
    }
};

struct SampleResultPoint {
    Vec3d position;
    float power_density = 0.0f;
};

struct MethodResult {
    std::string method;
    std::size_t target_points = 0;
    std::size_t valid_queries = 0;
    std::size_t projection_failures = 0;
    std::size_t visible_count = 0;
    std::size_t occluded_count = 0;
    std::size_t total_candidates = 0;
    float min_power_density = 0.0f;
    float max_power_density = 0.0f;
    double core_time_s = 0.0;
    double throughput_qps = 0.0;
    double avg_point_us = 0.0;
    double peak_memory_gb = 0.0;
};

double ParseDouble(const char* value, const std::string& name) {
    char* end = nullptr;
    const double parsed = std::strtod(value, &end);
    if (end == value || *end != '\0') {
        throw std::invalid_argument("Invalid numeric value for " + name + ": " + value);
    }
    return parsed;
}

int ParseInt(const char* value, const std::string& name) {
    char* end = nullptr;
    const long parsed = std::strtol(value, &end, 10);
    if (end == value || *end != '\0' || parsed < 0) {
        throw std::invalid_argument("Invalid integer value for " + name + ": " + value);
    }
    return static_cast<int>(parsed);
}

std::uint64_t ParseUInt64(const char* value, const std::string& name) {
    char* end = nullptr;
    const unsigned long long parsed = std::strtoull(value, &end, 10);
    if (end == value || *end != '\0') {
        throw std::invalid_argument("Invalid unsigned integer value for " + name + ": " + value);
    }
    return static_cast<std::uint64_t>(parsed);
}

std::size_t ParseSizeValue(const std::string& text) {
    char* end = nullptr;
    const unsigned long long parsed = std::strtoull(text.c_str(), &end, 10);
    if (end == text.c_str() || *end != '\0' || parsed == 0) {
        throw std::invalid_argument("Invalid point count: " + text);
    }
    return static_cast<std::size_t>(parsed);
}

std::vector<std::size_t> ParseSizes(const std::string& text) {
    std::vector<std::size_t> sizes;
    std::stringstream ss(text);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (!item.empty()) {
            sizes.push_back(ParseSizeValue(item));
        }
    }
    if (sizes.empty()) {
        throw std::invalid_argument("Point size list is empty");
    }
    return sizes;
}

void PrintUsage(const char* exe) {
    std::cout
        << "Usage: " << exe << " [options]\n\n"
        << "Methods:\n"
        << "  --method <vrpf|dsm|both>   Method to benchmark (default: both)\n"
        << "  --index <path>             RTree index file for VRPF\n"
        << "  --dem <path>               1 m DSM/DEM GeoTIFF for DSM method\n"
        << "  --sizes <a,b,c>            Target point counts (default: 250000,1000000,5000000,10000000)\n"
        << "  --threads <n>              OpenMP threads, 0 means runtime default\n"
        << "  --warmup-runs <n>          Untimed passes before measured run (default: 1)\n"
        << "  --csv <path>               CSV output path\n"
        << "  --no-csv                   Disable CSV output\n\n"
        << "Random 3D sampling:\n"
        << "  --min-lon/--max-lon <deg>  Longitude range\n"
        << "  --min-lat/--max-lat <deg>  Latitude range\n"
        << "  --min-alt/--max-alt <m>    Target absolute altitude range (default: 0..100)\n"
        << "  --seed <n>                 Random seed (default: 20260517)\n"
        << "  --fixed-height <m>         Compatibility alias: set min-alt=max-alt\n\n"
        << "Radar and RTree coordinates:\n"
        << "  --radar-lon <deg>          Radar longitude\n"
        << "  --radar-lat <deg>          Radar latitude\n"
        << "  --radar-alt <m>            Radar height above terrain\n"
        << "  --min-x <value>            EPSG:2326 min X used by index\n"
        << "  --min-y <value>            EPSG:2326 min Y used by index\n"
        << "  --index-range-x <value>    RTree X offset\n"
        << "  --index-range-y <value>    RTree Y offset\n";
}

bool ParseOptions(int argc, char** argv, Options& options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto value = [&](const std::string& name) -> const char* {
            if (i + 1 >= argc) {
                throw std::invalid_argument("Missing value for " + name);
            }
            return argv[++i];
        };

        if (arg == "--help" || arg == "-h") {
            options.show_help = true;
            return true;
        } else if (arg == "--method") {
            options.method = value(arg);
        } else if (arg == "--index") {
            options.index_file = value(arg);
        } else if (arg == "--dem") {
            options.dem_file = value(arg);
        } else if (arg == "--sizes") {
            options.sizes = ParseSizes(value(arg));
        } else if (arg == "--threads") {
            options.threads = ParseInt(value(arg), arg);
        } else if (arg == "--warmup-runs") {
            options.warmup_runs = ParseInt(value(arg), arg);
        } else if (arg == "--csv") {
            options.csv_file = value(arg);
        } else if (arg == "--no-csv") {
            options.write_csv = false;
        } else if (arg == "--min-lon") {
            options.bounds.min_lon = ParseDouble(value(arg), arg);
        } else if (arg == "--max-lon") {
            options.bounds.max_lon = ParseDouble(value(arg), arg);
        } else if (arg == "--min-lat") {
            options.bounds.min_lat = ParseDouble(value(arg), arg);
        } else if (arg == "--max-lat") {
            options.bounds.max_lat = ParseDouble(value(arg), arg);
        } else if (arg == "--fixed-height") {
            const double height = ParseDouble(value(arg), arg);
            options.min_alt_m = height;
            options.max_alt_m = height;
        } else if (arg == "--min-alt") {
            options.min_alt_m = ParseDouble(value(arg), arg);
        } else if (arg == "--max-alt") {
            options.max_alt_m = ParseDouble(value(arg), arg);
        } else if (arg == "--seed") {
            options.seed = ParseUInt64(value(arg), arg);
        } else if (arg == "--radar-lon") {
            options.radar_lon = ParseDouble(value(arg), arg);
        } else if (arg == "--radar-lat") {
            options.radar_lat = ParseDouble(value(arg), arg);
        } else if (arg == "--radar-alt") {
            options.radar_alt_relative_m = ParseDouble(value(arg), arg);
        } else if (arg == "--min-x") {
            options.coord.min_x = ParseDouble(value(arg), arg);
        } else if (arg == "--min-y") {
            options.coord.min_y = ParseDouble(value(arg), arg);
        } else if (arg == "--index-range-x") {
            options.coord.index_range_x = ParseDouble(value(arg), arg);
        } else if (arg == "--index-range-y") {
            options.coord.index_range_y = ParseDouble(value(arg), arg);
        } else {
            throw std::invalid_argument("Unknown option: " + arg);
        }
    }

    if (options.method != "vrpf" && options.method != "dsm" && options.method != "both") {
        throw std::invalid_argument("--method must be vrpf, dsm, or both");
    }
    if (options.bounds.min_lon >= options.bounds.max_lon ||
        options.bounds.min_lat >= options.bounds.max_lat) {
        throw std::invalid_argument("Invalid sampling bounds");
    }
    if (options.min_alt_m > options.max_alt_m) {
        throw std::invalid_argument("Invalid altitude range");
    }
    return true;
}

void ConfigureThreads(int threads) {
#ifdef _OPENMP
    if (threads > 0) {
        omp_set_num_threads(threads);
    }
#else
    (void)threads;
#endif
}

int ThreadCount() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

double CurrentPeakMemoryGB() {
    struct rusage usage {};
    if (getrusage(RUSAGE_SELF, &usage) != 0) {
        return 0.0;
    }
    return static_cast<double>(usage.ru_maxrss) / (1024.0 * 1024.0);
}

bool LonLatToRTreeCoord(const Projector& projector, double lon, double lat,
                        const RTreeCoordConfig& coord, double& rtree_x, double& rtree_y) {
    double x2326 = 0.0;
    double y2326 = 0.0;
    if (!projector.LonLatToProjected(lon, lat, x2326, y2326)) {
        return false;
    }
    rtree_x = coord.index_range_x + (x2326 - coord.min_x);
    rtree_y = coord.index_range_y + (y2326 - coord.min_y);
    return true;
}

std::uint64_t SplitMix64(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

double UnitRandomFromBits(std::uint64_t value) {
    return static_cast<double>(value >> 11) * (1.0 / 9007199254740992.0);
}

double Lerp(double a, double b, double t) {
    return a + (b - a) * t;
}

std::vector<SampleResultPoint> GenerateRandomSamples(const Options& options, std::size_t target_points) {
    std::vector<SampleResultPoint> results(target_points);
    const double lon_span_min = options.bounds.min_lon;
    const double lon_span_max = options.bounds.max_lon;
    const double lat_span_min = options.bounds.min_lat;
    const double lat_span_max = options.bounds.max_lat;
    const double alt_span_min = options.min_alt_m;
    const double alt_span_max = options.max_alt_m;

    #pragma omp parallel for schedule(static)
    for (std::size_t idx = 0; idx < target_points; ++idx) {
        const std::uint64_t base = options.seed ^ (static_cast<std::uint64_t>(idx) * 0x9e3779b97f4a7c15ULL);
        const double lon_u = UnitRandomFromBits(SplitMix64(base));
        const double lat_u = UnitRandomFromBits(SplitMix64(base + 0xbf58476d1ce4e5b9ULL));
        const double alt_u = UnitRandomFromBits(SplitMix64(base + 0x94d049bb133111ebULL));
        results[idx].position.x = Lerp(lon_span_min, lon_span_max, lon_u);
        results[idx].position.y = Lerp(lat_span_min, lat_span_max, lat_u);
        results[idx].position.z = Lerp(alt_span_min, alt_span_max, alt_u);
        results[idx].power_density = 0.0f;
    }
    return results;
}

inline bool LineSegmentIntersectsTriangleFast(const ValueType& id, const double a[3], const double b[3]) {
    const double v0x = std::get<0>(id);
    const double v0y = std::get<1>(id);
    const double v0z = std::get<2>(id);
    const double v1x = std::get<3>(id);
    const double v1y = std::get<4>(id);
    const double v1z = std::get<5>(id);
    const double v2x = std::get<6>(id);
    const double v2y = std::get<7>(id);
    const double v2z = std::get<8>(id);
    const double dirx = b[0] - a[0];
    const double diry = b[1] - a[1];
    const double dirz = b[2] - a[2];
    const double e1x = v1x - v0x;
    const double e1y = v1y - v0y;
    const double e1z = v1z - v0z;
    const double e2x = v2x - v0x;
    const double e2y = v2y - v0y;
    const double e2z = v2z - v0z;
    const double px = diry * e2z - dirz * e2y;
    const double py = dirz * e2x - dirx * e2z;
    const double pz = dirx * e2y - diry * e2x;
    const double det = e1x * px + e1y * py + e1z * pz;
    if (det > -EPSILON && det < EPSILON) {
        return false;
    }
    const double inv_det = 1.0 / det;
    const double tx = a[0] - v0x;
    const double ty = a[1] - v0y;
    const double tz = a[2] - v0z;
    const double u = (tx * px + ty * py + tz * pz) * inv_det;
    if (u < 0.0 || u > 1.0) {
        return false;
    }
    const double qx = ty * e1z - tz * e1y;
    const double qy = tz * e1x - tx * e1z;
    const double qz = tx * e1y - ty * e1x;
    const double v = (dirx * qx + diry * qy + dirz * qz) * inv_det;
    if (v < 0.0 || u + v > 1.0) {
        return false;
    }
    const double t = (e2x * qx + e2y * qy + e2z * qz) * inv_det;
    return t > EPSILON && t < 1.0 - EPSILON;
}

float CalculatePowerDensityForTarget(Radar& radar, const Vec3d& target, caltools& ct) {
    return radar.CalculateSinglePointPowerDensity(target, &ct);
}

MethodResult RunVrpf(Radar& radar, RTree3d& rtree, std::vector<SampleResultPoint>& results, std::size_t count,
                     const double radar_point[3], const RTreeCoordConfig& coord, int warmup_runs) {
    auto run_once = [&](bool measured) {
        std::size_t total_candidates = 0;
        std::size_t occluded_count = 0;
        std::size_t projection_failures = 0;
        float min_density = std::numeric_limits<float>::max();
        float max_density = 0.0f;
        const auto start = std::chrono::steady_clock::now();

        #pragma omp parallel reduction(+:total_candidates, occluded_count, projection_failures) \
                             reduction(min:min_density) reduction(max:max_density)
        {
            Projector projector;
            const bool projector_ok = projector.IsValid();
            caltools ct;
            std::size_t local_candidates = 0;
            #pragma omp for schedule(guided, 512)
            for (std::size_t i = 0; i < count; ++i) {
                SampleResultPoint& point = results[i];
                point.power_density = 0.0f;
                if (!projector_ok) {
                    ++projection_failures;
                    continue;
                }
                double rtree_x = 0.0;
                double rtree_y = 0.0;
                if (!LonLatToRTreeCoord(projector, point.position.x, point.position.y, coord, rtree_x, rtree_y)) {
                    ++projection_failures;
                    continue;
                }
                local_candidates = 0;
                const double target[3] = {rtree_x, rtree_y, point.position.z};
                auto callback = [&local_candidates](const ValueType& id, const double a[3], const double b[3]) -> bool {
                    ++local_candidates;
                    return !LineSegmentIntersectsTriangleFast(id, a, b);
                };
                const bool occluded = rtree.Intersect3dFast(radar_point, target, callback);
                total_candidates += local_candidates;
                if (occluded) {
                    ++occluded_count;
                    continue;
                }
                const float density = CalculatePowerDensityForTarget(radar, point.position, ct);
                point.power_density = density;
                if (density > 0.0f) {
                    min_density = std::min(min_density, density);
                    max_density = std::max(max_density, density);
                }
            }
        }

        const auto end = std::chrono::steady_clock::now();
        MethodResult result;
        result.method = "VRPF";
        result.target_points = count;
        result.projection_failures = projection_failures;
        result.valid_queries = count - projection_failures;
        result.total_candidates = total_candidates;
        result.occluded_count = occluded_count;
        result.visible_count = result.valid_queries - occluded_count;
        result.min_power_density = min_density == std::numeric_limits<float>::max() ? 0.0f : min_density;
        result.max_power_density = max_density;
        result.core_time_s = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
        if (result.valid_queries != 0 && result.core_time_s > 0.0) {
            result.throughput_qps = static_cast<double>(result.valid_queries) / result.core_time_s;
            result.avg_point_us = result.core_time_s * 1000000.0 / static_cast<double>(result.valid_queries);
        }
        result.peak_memory_gb = CurrentPeakMemoryGB();
        (void)measured;
        return result;
    };

    for (int i = 0; i < warmup_runs; ++i) {
        MethodResult warmup = run_once(false);
        std::cout << "  VRPF warmup " << (i + 1) << ": " << warmup.core_time_s
                  << " s, " << warmup.throughput_qps << " q/s" << std::endl;
    }
    return run_once(true);
}

float DemHeightAt(const hiradar::DEMLoader& dem, double lon, double lat) {
    return dem.GetHeight(lon, lat);
}

bool IsOccludedDsmFast(const hiradar::DEMLoader& dem, const Vec3d& radar, const Vec3d& target) {
    const double dx_deg = target.x - radar.x;
    const double dy_deg = target.y - radar.y;
    if (std::abs(dx_deg) < 1e-12 && std::abs(dy_deg) < 1e-12) {
        return false;
    }

    const double center_lat = 0.5 * (radar.y + target.y);
    const double lat_rad = center_lat * M_PI / 180.0;
    const double meters_per_deg_lat = 111320.0;
    const double meters_per_deg_lon = meters_per_deg_lat * std::cos(lat_rad);
    const double dx_m = dx_deg * meters_per_deg_lon;
    const double dy_m = dy_deg * meters_per_deg_lat;
    const double horizontal_dist_m = std::sqrt(dx_m * dx_m + dy_m * dy_m);
    const double dem_resolution_m = std::abs(dem.GetResolution()) * meters_per_deg_lat;
    int sample_count = static_cast<int>(std::ceil(horizontal_dist_m / dem_resolution_m));
    sample_count = std::max(sample_count, 2);

    constexpr float terrain_clearance = 0.05f;
    const float nodata = dem.GetNoDataValue();
    for (int i = 1; i < sample_count; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(sample_count);
        const double lon = radar.x + t * dx_deg;
        const double lat = radar.y + t * dy_deg;
        const double ray_h = radar.z + t * (target.z - radar.z);
        const float terrain_h = DemHeightAt(dem, lon, lat);
        if (terrain_h == nodata) {
            continue;
        }
        if (ray_h < static_cast<double>(terrain_h - terrain_clearance)) {
            return true;
        }
    }
    return false;
}

MethodResult RunDsm(Radar& radar_model, const hiradar::DEMLoader& dem, std::vector<SampleResultPoint>& results,
                    std::size_t count, const Vec3d& radar, int warmup_runs) {
    auto run_once = [&]() {
        std::size_t occluded_count = 0;
        float min_density = std::numeric_limits<float>::max();
        float max_density = 0.0f;
        const auto start = std::chrono::steady_clock::now();

        #pragma omp parallel reduction(+:occluded_count) reduction(min:min_density) reduction(max:max_density)
        {
            caltools ct;
            #pragma omp for schedule(guided, 512)
            for (std::size_t i = 0; i < count; ++i) {
                SampleResultPoint& point = results[i];
                point.power_density = 0.0f;
                if (IsOccludedDsmFast(dem, radar, point.position)) {
                    ++occluded_count;
                    continue;
                }
                const float density = CalculatePowerDensityForTarget(radar_model, point.position, ct);
                point.power_density = density;
                if (density > 0.0f) {
                    min_density = std::min(min_density, density);
                    max_density = std::max(max_density, density);
                }
            }
        }

        const auto end = std::chrono::steady_clock::now();
        MethodResult result;
        result.method = "1m_DSM";
        result.target_points = count;
        result.valid_queries = count;
        result.occluded_count = occluded_count;
        result.visible_count = count - occluded_count;
        result.min_power_density = min_density == std::numeric_limits<float>::max() ? 0.0f : min_density;
        result.max_power_density = max_density;
        result.core_time_s = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
        result.throughput_qps = static_cast<double>(count) / result.core_time_s;
        result.avg_point_us = result.core_time_s * 1000000.0 / static_cast<double>(count);
        result.peak_memory_gb = CurrentPeakMemoryGB();
        return result;
    };

    for (int i = 0; i < warmup_runs; ++i) {
        MethodResult warmup = run_once();
        std::cout << "  DSM warmup " << (i + 1) << ": " << warmup.core_time_s
                  << " s, " << warmup.throughput_qps << " q/s" << std::endl;
    }
    return run_once();
}

void PrintResult(const MethodResult& result) {
    std::cout << std::fixed << std::setprecision(6)
              << result.method
              << " points=" << result.target_points
              << " time_s=" << result.core_time_s
              << " qps=" << result.throughput_qps
              << " avg_us=" << result.avg_point_us
              << " peak_gb=" << result.peak_memory_gb
              << " visible=" << result.visible_count
              << " occluded=" << result.occluded_count
              << " min_density=" << result.min_power_density
              << " max_density=" << result.max_power_density;
    if (result.method == "VRPF") {
        std::cout << " candidates=" << result.total_candidates;
    }
    std::cout << std::endl;
}

void WriteCsv(const std::vector<MethodResult>& results, const std::string& path) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Cannot open CSV: " + path);
    }
    out << "method,target_points,valid_queries,projection_failures,core_time_s,throughput_qps,avg_point_us,peak_memory_gb,"
           "visible_count,occluded_count,total_candidates,min_power_density,max_power_density\n";
    out << std::fixed << std::setprecision(6);
    for (const MethodResult& r : results) {
        out << r.method << ','
            << r.target_points << ','
            << r.valid_queries << ','
            << r.projection_failures << ','
            << r.core_time_s << ','
            << r.throughput_qps << ','
            << r.avg_point_us << ','
            << r.peak_memory_gb << ','
            << r.visible_count << ','
            << r.occluded_count << ','
            << r.total_candidates << ','
            << r.min_power_density << ','
            << r.max_power_density << '\n';
    }
}

} // namespace

int main(int argc, char** argv) {
    Options options;
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

    ConfigureThreads(options.threads);
    std::sort(options.sizes.begin(), options.sizes.end());
    const std::size_t max_count = options.sizes.back();

    std::cout << "========================================\n";
    std::cout << "  VRPF vs 1 m DSM Scaling Benchmark\n";
    std::cout << "========================================\n";
    std::cout << "Method: " << options.method << "\n";
    std::cout << "Threads: " << ThreadCount() << "\n";
    std::cout << "Max target points: " << max_count << "\n";
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Lon range: [" << options.bounds.min_lon << ", " << options.bounds.max_lon << "]\n";
    std::cout << "Lat range: [" << options.bounds.min_lat << ", " << options.bounds.max_lat << "]\n";
    std::cout << "Alt range: [" << options.min_alt_m << ", " << options.max_alt_m << "] m\n";
    std::cout << "Random seed: " << options.seed << "\n";
    std::cout << std::defaultfloat;

    Projector projector;
    if (!projector.IsValid()) {
        std::cerr << "Error: failed to create EPSG:4326 -> EPSG:2326 projector\n";
        return 1;
    }

    double radar_rtree_x = 0.0;
    double radar_rtree_y = 0.0;
    if (!LonLatToRTreeCoord(projector, options.radar_lon, options.radar_lat, options.coord,
                            radar_rtree_x, radar_rtree_y)) {
        std::cerr << "Error: failed to project radar position\n";
        return 1;
    }

    Position radar_pos_dms;
    DecimalToDMS(radar_pos_dms.lon, static_cast<float>(options.radar_lon));
    DecimalToDMS(radar_pos_dms.lat, static_cast<float>(options.radar_lat));
    radar_pos_dms.alt = static_cast<float>(options.radar_alt_relative_m);
    const float base_phi = 0.0f;
    const float base_theta = 0.0f;
    std::unique_ptr<Radar> radar(CityGuardRadar(radar_pos_dms, base_phi, base_theta));
    radar->_cached_abs_alt_m = options.radar_alt_relative_m;

    std::vector<MethodResult> results;

    RTree3d rtree;
    double radar_abs_alt_m = options.radar_alt_relative_m;
    double radar_point[3] = {radar_rtree_x, radar_rtree_y, radar_abs_alt_m};
    if (options.method == "vrpf" || options.method == "both") {
        std::cout << "Loading VRPF RTree: " << options.index_file << std::endl;
        if (!rtree.Load(options.index_file.c_str())) {
            std::cerr << "Error: cannot load RTree index\n";
            return 1;
        }
        radar->BindRtree(&rtree);
        radar->Move(radar_pos_dms, base_phi, base_theta, 0.0f);
        const double terrain_h = rtree.Getheight3d(radar_rtree_x, radar_rtree_y, HeightCallback);
        if (terrain_h >= MINH) {
            radar_abs_alt_m = terrain_h + options.radar_alt_relative_m;
            radar_point[2] = radar_abs_alt_m;
        }
        radar->_cached_abs_alt_m = radar_abs_alt_m;
        std::cout << "VRPF radar RTree xyz: " << radar_point[0] << ", "
                  << radar_point[1] << ", " << radar_point[2] << std::endl;
        std::cout << "Peak memory after RTree load: " << CurrentPeakMemoryGB() << " GB\n";
    }

    hiradar::DEMLoader dem;
    Vec3d radar_lonlat;
    radar_lonlat.x = options.radar_lon;
    radar_lonlat.y = options.radar_lat;
    radar_lonlat.z = radar_abs_alt_m;
    if (options.method == "dsm" || options.method == "both") {
        std::cout << "Loading 1 m DSM: " << options.dem_file << std::endl;
        if (!dem.Load(options.dem_file)) {
            std::cerr << "Error: cannot load DSM file\n";
            return 1;
        }
        if (options.method == "dsm") {
            const float dem_radar_h = dem.GetHeight(options.radar_lon, options.radar_lat);
            if (dem_radar_h != dem.GetNoDataValue()) {
                radar_lonlat.z = dem_radar_h + options.radar_alt_relative_m;
                radar->_cached_abs_alt_m = radar_lonlat.z;
            }
        }
        std::cout << "DSM radar lon/lat/z: " << radar_lonlat.x << ", "
                  << radar_lonlat.y << ", " << radar_lonlat.z << std::endl;
        std::cout << "Peak memory after DSM load: " << CurrentPeakMemoryGB() << " GB\n";
    }

    for (std::size_t count : options.sizes) {
        std::cout << "\n--- target points: " << count << " ---\n";
        const auto sample_start = std::chrono::steady_clock::now();
        std::vector<SampleResultPoint> sample_results = GenerateRandomSamples(options, count);
        const auto sample_end = std::chrono::steady_clock::now();
        const double sample_time_s =
            std::chrono::duration_cast<std::chrono::duration<double>>(sample_end - sample_start).count();
        std::cout << "Random 3D sample points generated: " << sample_results.size()
                  << " in " << sample_time_s << " s"
                  << ". Peak memory now: " << CurrentPeakMemoryGB() << " GB\n";
        if (options.method == "vrpf" || options.method == "both") {
            MethodResult vrpf = RunVrpf(*radar, rtree, sample_results, count, radar_point, options.coord, options.warmup_runs);
            PrintResult(vrpf);
            results.push_back(vrpf);
        }
        if (options.method == "dsm" || options.method == "both") {
            MethodResult dsm = RunDsm(*radar, dem, sample_results, count, radar_lonlat, options.warmup_runs);
            PrintResult(dsm);
            results.push_back(dsm);
        }
    }

    if (options.method == "both") {
        std::cout << "\n========== Table Rows ==========\n";
        for (std::size_t count : options.sizes) {
            const MethodResult* vrpf = nullptr;
            const MethodResult* dsm = nullptr;
            for (const MethodResult& result : results) {
                if (result.target_points == count && result.method == "VRPF") {
                    vrpf = &result;
                } else if (result.target_points == count && result.method == "1m_DSM") {
                    dsm = &result;
                }
            }
            if (vrpf && dsm) {
                const double speedup = dsm->core_time_s / vrpf->core_time_s;
                std::cout << "points=" << count
                          << " vrpf_s=" << vrpf->core_time_s
                          << " dsm_s=" << dsm->core_time_s
                          << " speedup=" << speedup
                          << " vrpf_peak_gb=" << vrpf->peak_memory_gb
                          << " dsm_peak_gb=" << dsm->peak_memory_gb
                          << '\n';
            }
        }
    }

    if (options.write_csv) {
        WriteCsv(results, options.csv_file);
        std::cout << "\nCSV written: " << options.csv_file << std::endl;
    }

    std::cout << "\nFinal process peak memory: " << CurrentPeakMemoryGB() << " GB\n";
    return 0;
}
