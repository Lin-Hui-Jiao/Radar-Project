#include "hiradar/dem_generator.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace hiradar {

void DEMGenerator::SetGeoTransformEPSG2326(
    GDALDataset* dataset,
    const RTreeCoordConfig& coord_config,
    double resolution_m
) {
    const double top_y = coord_config.min_y + coord_config.area_height_m;
    double geotransform[6] = {
        coord_config.min_x,
        resolution_m,
        0.0,
        top_y,
        0.0,
        -resolution_m
    };
    dataset->SetGeoTransform(geotransform);

    OGRSpatialReference srs;
    srs.importFromEPSG(2326);
    char* wkt = nullptr;
    srs.exportToWkt(&wkt);
    dataset->SetProjection(wkt);
    CPLFree(wkt);
}

float DEMGenerator::QueryHeightAtProjectedCoord(
    RTree3d* rtree,
    double x2326,
    double y2326,
    const RTreeCoordConfig& coord_config
) {
    if (!rtree) {
        return -9999.0f;
    }

    const double rtree_x = coord_config.index_range_x + (x2326 - coord_config.min_x);
    const double rtree_y = coord_config.index_range_y + (y2326 - coord_config.min_y);
    const double height = rtree->Getheight3d(rtree_x, rtree_y, HeightCallback);

    if (height < MINH) {
        return -9999.0f;
    }
    return static_cast<float>(height);
}

bool DEMGenerator::GenerateDEMFromTriangles(
    RTree3d* rtree,
    double lon_range[2],
    double lat_range[2],
    double resolution_m,
    const char* output_path,
    const RTreeCoordConfig& coord_config
) {
    (void)lon_range;
    (void)lat_range;

    if (!rtree) {
        std::cerr << "Error: RTree is null" << std::endl;
        return false;
    }
    if (resolution_m <= 0.0) {
        std::cerr << "Error: invalid DEM resolution: " << resolution_m << std::endl;
        return false;
    }
    if (coord_config.area_width_m <= 0.0 || coord_config.area_height_m <= 0.0) {
        std::cerr << "Error: invalid projected area size: "
                  << coord_config.area_width_m << " x "
                  << coord_config.area_height_m << std::endl;
        return false;
    }

    const int width = static_cast<int>(std::ceil(coord_config.area_width_m / resolution_m));
    const int height = static_cast<int>(std::ceil(coord_config.area_height_m / resolution_m));

    if (width <= 0 || height <= 0) {
        std::cerr << "Error: invalid DEM raster size: " << width << " x " << height << std::endl;
        return false;
    }

    std::cout << "\n========================================" << std::endl;
    std::cout << "Generating DEM from triangulated surface" << std::endl;
    std::cout << "Output CRS: EPSG:2326 (Hong Kong 1980 Grid System)" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Output: " << output_path << std::endl;
    std::cout << "Resolution: " << resolution_m << " m/pixel" << std::endl;
    std::cout << "Projected X range: [" << coord_config.min_x << ", "
              << (coord_config.min_x + coord_config.area_width_m) << "] m" << std::endl;
    std::cout << "Projected Y range: [" << coord_config.min_y << ", "
              << (coord_config.min_y + coord_config.area_height_m) << "] m" << std::endl;
    std::cout << "Raster size: " << width << " x " << height << " pixels" << std::endl;

    GDALAllRegister();

    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GTiff");
    if (!driver) {
        std::cerr << "Error: GTiff driver not available" << std::endl;
        return false;
    }

    GDALDataset* dataset = driver->Create(output_path, width, height, 1, GDT_Float32, nullptr);
    if (!dataset) {
        std::cerr << "Error: Failed to create output dataset" << std::endl;
        return false;
    }

    SetGeoTransformEPSG2326(dataset, coord_config, resolution_m);

    GDALRasterBand* band = dataset->GetRasterBand(1);
    band->SetNoDataValue(-9999.0f);

    std::cout << "Querying heights from R-Tree (EPSG:2326)..." << std::endl;

    float* scanline = new float[width];
    const int progress_interval = std::max(1, height / 20);
    const double top_y = coord_config.min_y + coord_config.area_height_m;

    size_t nodata_count = 0;
    for (int row = 0; row < height; row++) {
        if (row % progress_interval == 0) {
            std::cout << "Progress: " << (row * 100 / height) << "%" << std::endl;
        }

        const double y2326 = top_y - (static_cast<double>(row) + 0.5) * resolution_m;
        for (int col = 0; col < width; col++) {
            const double x2326 = coord_config.min_x + (static_cast<double>(col) + 0.5) * resolution_m;
            scanline[col] = QueryHeightAtProjectedCoord(rtree, x2326, y2326, coord_config);
            if (scanline[col] == -9999.0f) {
                nodata_count++;
            }
        }

        CPLErr err = band->RasterIO(
            GF_Write, 0, row, width, 1, scanline, width, 1, GDT_Float32, 0, 0
        );
        if (err != CE_None) {
            std::cerr << "Error: Failed to write raster data at row " << row << std::endl;
        }
    }

    delete[] scanline;
    GDALClose(dataset);

    std::cout << "Progress: 100%" << std::endl;
    std::cout << "NoData pixels: " << nodata_count << " / "
              << static_cast<size_t>(width) * static_cast<size_t>(height) << std::endl;
    std::cout << "DEM generation completed: " << output_path << std::endl;
    std::cout << "========================================\n" << std::endl;

    return true;
}

void DEMGenerator::GenerateMultiResolutionDEMs(
    RTree3d* rtree,
    double lon_range[2],
    double lat_range[2],
    const double* resolutions_m,
    int num_resolutions,
    const char* output_prefix,
    const RTreeCoordConfig& coord_config
) {
    std::cout << "\n****************************************" << std::endl;
    std::cout << "Generating Multi-Resolution DEM Set" << std::endl;
    std::cout << "****************************************" << std::endl;
    std::cout << "Number of resolutions: " << num_resolutions << std::endl;

    for (int i = 0; i < num_resolutions; i++) {
        const double res_m = resolutions_m[i];

        std::string filename = std::string(output_prefix);
        if (res_m < 1.0) {
            filename += std::to_string(res_m).substr(0, 3) + "m.tif";
        } else {
            filename += std::to_string(static_cast<int>(res_m)) + "m.tif";
        }

        std::cout << "\n[" << (i + 1) << "/" << num_resolutions << "] Generating "
                  << filename << " (resolution: " << res_m << " m)" << std::endl;

        bool success = GenerateDEMFromTriangles(
            rtree, lon_range, lat_range, res_m, filename.c_str(), coord_config
        );

        if (!success) {
            std::cerr << "Warning: Failed to generate " << filename << std::endl;
        }
    }

    std::cout << "\n****************************************" << std::endl;
    std::cout << "Multi-Resolution DEM Set Completed!" << std::endl;
    std::cout << "****************************************\n" << std::endl;
}

} // namespace hiradar
