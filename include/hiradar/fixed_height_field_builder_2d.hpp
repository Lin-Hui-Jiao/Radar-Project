#pragma once

#include <string>

#include "hiradar/energy_grid_2d.hpp"

namespace hiradar {

struct MeshFieldBuildConfig {
    double lon_min = 114.164571;
    double lon_max = 114.169422;
    double lat_min = 22.278365;
    double lat_max = 22.28288;
    double fixed_height_m = 50.0;
    double resolution_m = 1.0;
    double radar_lon_deg = 114.1670;
    double radar_lat_deg = 22.2806;
    double radar_relative_alt_m = 80.0;
    float base_phi = 0.0f;
    float base_theta = 0.0f;
    std::string rtree_file = "../../test_area.3idx";
};

class FixedHeightFieldBuilder2D {
public:
    static EnergyGrid2D BuildFromMesh(const MeshFieldBuildConfig& config);
};

}  // namespace hiradar
