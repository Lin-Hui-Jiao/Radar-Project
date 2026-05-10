#pragma once

#include "hiradar/RTree.h"
#include "hiradar/dataset_config.hpp"
#include "hiradar/dem_loader.hpp"
#include "hiradar/dem_visibility.hpp"
#include "hiradar/occlusion_utils.h"

#include <string>
#include <vector>

class Radar;

namespace hiradar {

class FixedHeightComparator {
public:
    struct GridPoint {
        Vec3d position;      // EPSG:4326: x=lon, y=lat, z=alt ASL.
        double x2326 = 0.0;  // EPSG:2326 easting, meters.
        double y2326 = 0.0;  // EPSG:2326 northing, meters.
        double rtree_x = 0.0;
        double rtree_y = 0.0;
        bool visible_mesh = false;
        bool visible_dem = false;
        float power_density = 0.0f;
    };

    struct ExperimentResult {
        std::vector<GridPoint> grid_points;

        size_t total_points = 0;
        size_t visible_mesh_count = 0;
        size_t visible_dem_count = 0;
        size_t disagreement_count = 0;

        int grid_width = 0;
        int grid_height = 0;

        float min_power_density = 0.0f;
        float max_power_density = 0.0f;

        double fixed_height = 0.0;
        double grid_resolution = 0.0;
        std::string dem_file;
        double dem_resolution = 0.0;

        double total_time = 0.0;
        double mesh_visibility_time = 0.0;
        double dem_visibility_time = 0.0;
        double power_calculation_time = 0.0;
        int num_threads = 0;
    };

    static ExperimentResult RunExperiment(
        Radar* radar,
        const Vec3d& radar_pos,
        double lon_range[2],
        double lat_range[2],
        double fixed_height,
        double grid_resolution_m,
        RTree3d* rtree,
        DEMLoader* dem,
        const RTreeCoordConfig& coord_config
    );

    static void GenerateVisualization(
        const ExperimentResult& result,
        const std::string& output_file
    );

    static void ExportResultsToCSV(
        const ExperimentResult& result,
        const std::string& output_file
    );

    static void PrintSummary(const ExperimentResult& result);
};

} // namespace hiradar
