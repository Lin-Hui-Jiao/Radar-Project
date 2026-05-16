#pragma once

namespace hiradar {

struct DatasetBounds {
    double min_lon = 114.174274;
    double min_lat = 22.262560;
    double max_lon = 114.193680;
    double max_lat = 22.280623;
};

struct RTreeCoordConfig {
    double index_range_x = -974.6354;
    double index_range_y = -1706.5722;
    double min_x = 836000.0;
    double min_y = 813750.0;
    double area_width_m = 2000.0;
    double area_height_m = 2000.0;
};

} // namespace hiradar



// lon = 114.1838991865
// lat = 22.2738043465






// ./fixed_height_experiment --index ../../test_area_2_2.3idx --min-x 836000 --min-y 813750 --index-range-x -974.6354 --index-range-y -1706.5722 
// --radar-lon 114.183899 --radar-lat 22.273804 --radar-alt 80 --min-lon 114.174274 --max-lon 114.193680 --min-lat 22.262560 --max-lat 22.280623 --fixed-height 30