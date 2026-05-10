#pragma once

namespace hiradar {

struct DatasetBounds {
    double min_lon = 114.157292795;
    double min_lat = 22.273848477;
    double max_lon = 114.166996872;
    double max_lat = 22.282879994;
};

struct RTreeCoordConfig {
    double index_range_x = -2724.635391;
    double index_range_y = -456.572217;
    double min_x = 834250.0;
    double min_y = 815000.0;
    double area_width_m = 1000.0;
    double area_height_m = 1000.0;
};

} // namespace hiradar
