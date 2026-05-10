#pragma once

#include "hiradar/RTree.h"
#include "hiradar/dataset_config.hpp"
#include "hiradar/occlusion_utils.h"

#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include <string>

namespace hiradar {

class DEMGenerator {
public:
    static bool GenerateDEMFromTriangles(
        RTree3d* rtree,
        double lon_range[2],
        double lat_range[2],
        double resolution_m,
        const char* output_path,
        const RTreeCoordConfig& coord_config = RTreeCoordConfig{}
    );

    static void GenerateMultiResolutionDEMs(
        RTree3d* rtree,
        double lon_range[2],
        double lat_range[2],
        const double* resolutions_m,
        int num_resolutions,
        const char* output_prefix,
        const RTreeCoordConfig& coord_config = RTreeCoordConfig{}
    );

private:
    static void SetGeoTransformEPSG2326(
        GDALDataset* dataset,
        const RTreeCoordConfig& coord_config,
        double resolution_m
    );

    static float QueryHeightAtProjectedCoord(
        RTree3d* rtree,
        double x2326,
        double y2326,
        const RTreeCoordConfig& coord_config
    );
};

} // namespace hiradar
