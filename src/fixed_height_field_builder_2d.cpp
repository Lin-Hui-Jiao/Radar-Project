#include "hiradar/fixed_height_field_builder_2d.hpp"

#include <algorithm>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <vector>

#include <proj.h>

#include "hiradar/grid.hpp"
#include "hiradar/radar.hpp"
#include "hiradar/radar_pool.hpp"

namespace hiradar {

namespace {

constexpr double kIndexMinX = 835000.0;
constexpr double kIndexMinY = 815500.0;
constexpr double kIndexRangeX = -1975.0;
constexpr double kIndexRangeY = 43.0;

class IndexProjection {
public:
    IndexProjection() {
        context_ = proj_context_create();
        transform_ = proj_create_crs_to_crs(context_, "EPSG:4326", "EPSG:2326", nullptr);
        if (transform_ == nullptr) {
            throw std::runtime_error("Failed to create EPSG:4326 -> EPSG:2326 projection");
        }
    }

    ~IndexProjection() {
        if (transform_ != nullptr) {
            proj_destroy(transform_);
        }
        if (context_ != nullptr) {
            proj_context_destroy(context_);
        }
    }

    Vec3d ToIndexLocal(double lon_deg, double lat_deg) const {
        const PJ_COORD coord = proj_coord(lat_deg, lon_deg, 0.0, 0.0);
        const PJ_COORD projected = proj_trans(transform_, PJ_FWD, coord);
        return Vec3d(projected.xy.y - kIndexMinX, projected.xy.x - kIndexMinY, 0.0);
    }

private:
    PJ_CONTEXT* context_ = nullptr;
    PJ* transform_ = nullptr;
};

std::string ResolveExistingFile(const std::string& configured_path) {
    const std::vector<std::string> candidates = {
        configured_path,
        "../../test_area.3idx",
        "../test_area.3idx",
        "test_area.3idx",
    };

    for (const std::string& candidate : candidates) {
        if (!candidate.empty() && std::filesystem::exists(candidate)) {
            return candidate;
        }
    }

    throw std::runtime_error("Cannot locate RTree index file. Checked configured path and common defaults");
}

bool IsOccludedByMesh(
    const Vec3d& radar_local_xy,
    double radar_abs_alt_m,
    const Vec3d& target_local_xy,
    double target_abs_alt_m,
    const RTree3d& rtree
) {
    const double start[3] = {
        kIndexRangeX + radar_local_xy.x,
        kIndexRangeY + radar_local_xy.y,
        radar_abs_alt_m,
    };
    const double end[3] = {
        kIndexRangeX + target_local_xy.x,
        kIndexRangeY + target_local_xy.y,
        target_abs_alt_m,
    };
    return rtree.Intersect3d(start, end, IntersectCallback);
}

}  // namespace

EnergyGrid2D FixedHeightFieldBuilder2D::BuildFromMesh(const MeshFieldBuildConfig& config) {
    if (config.lon_min >= config.lon_max || config.lat_min >= config.lat_max) {
        throw std::runtime_error("Mesh field builder requires valid lon/lat ranges");
    }
    if (config.resolution_m <= 0.0) {
        throw std::runtime_error("Mesh field builder requires a positive grid resolution");
    }

    RTree3d rtree;
    const std::string rtree_path = ResolveExistingFile(config.rtree_file);
    if (!rtree.Load(rtree_path.c_str())) {
        throw std::runtime_error("Failed to load RTree index file: " + rtree_path);
    }

    Position radar_position{};
    DecimalToDMS(radar_position.lon, static_cast<float>(config.radar_lon_deg));
    DecimalToDMS(radar_position.lat, static_cast<float>(config.radar_lat_deg));
    radar_position.alt = static_cast<float>(config.radar_relative_alt_m);

    std::unique_ptr<Radar> radar(CityGuardRadar(radar_position, config.base_phi, config.base_theta));
    radar->BindRtree(&rtree);
    radar->Move(radar_position, config.base_phi, config.base_theta, 0.0f);

    const double center_lon = 0.5 * (config.lon_min + config.lon_max);
    const double center_lat = 0.5 * (config.lat_min + config.lat_max);
    LocalTangentFrame frame(center_lon, center_lat, config.fixed_height_m);

    const Vec3d corners[4] = {
        frame.ToENU(config.lon_min, config.lat_min, config.fixed_height_m),
        frame.ToENU(config.lon_max, config.lat_min, config.fixed_height_m),
        frame.ToENU(config.lon_min, config.lat_max, config.fixed_height_m),
        frame.ToENU(config.lon_max, config.lat_max, config.fixed_height_m),
    };

    double x_min_m = corners[0].x;
    double x_max_m = corners[0].x;
    double y_min_m = corners[0].y;
    double y_max_m = corners[0].y;
    for (const Vec3d& corner : corners) {
        x_min_m = std::min(x_min_m, corner.x);
        x_max_m = std::max(x_max_m, corner.x);
        y_min_m = std::min(y_min_m, corner.y);
        y_max_m = std::max(y_max_m, corner.y);
    }

    const int width = static_cast<int>(std::floor((x_max_m - x_min_m) / config.resolution_m + 0.5)) + 1;
    const int height = static_cast<int>(std::floor((y_max_m - y_min_m) / config.resolution_m + 0.5)) + 1;
    if (width <= 0 || height <= 0) {
        throw std::runtime_error("Mesh field builder recovered invalid grid dimensions");
    }

    EnergyGrid2D grid;
    grid.width = width;
    grid.height = height;
    grid.resolution_m = config.resolution_m;
    grid.fixed_height_m = config.fixed_height_m;
    grid.origin_frame = frame;
    grid.x_min_m = x_min_m;
    grid.y_min_m = y_min_m;
    grid.values.assign(static_cast<size_t>(width) * static_cast<size_t>(height), 0.0f);

    IndexProjection projection;
    const Vec3d radar_local_xy = projection.ToIndexLocal(config.radar_lon_deg, config.radar_lat_deg);
    caltools calculator;

    for (int iy = 0; iy < height; ++iy) {
        const double north_m = y_min_m + static_cast<double>(iy) * config.resolution_m;
        for (int ix = 0; ix < width; ++ix) {
            const double east_m = x_min_m + static_cast<double>(ix) * config.resolution_m;
            const Vec3d geodetic = frame.FromENU(east_m, north_m, 0.0);
            const Vec3d target_geo(geodetic.x, geodetic.y, config.fixed_height_m);
            const Vec3d target_local_xy = projection.ToIndexLocal(target_geo.x, target_geo.y);

            const bool occluded = IsOccludedByMesh(
                radar_local_xy,
                radar->_cached_abs_alt_m,
                target_local_xy,
                config.fixed_height_m,
                rtree
            );

            if (!occluded) {
                grid.values[static_cast<size_t>(iy) * width + ix] =
                    radar->CalculateSinglePointPowerDensity(target_geo, &calculator);
            }
        }
    }

    return grid;
}

}  // namespace hiradar
