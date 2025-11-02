#include "hiradar/dem_generator.hpp"
#include "hiradar/caltools.h"
#include <proj.h>
#include <iostream>
#include <cmath>
#include <limits>

namespace hiradar {

// 坐标转换常量（参考radar.hpp line 52-76）
const double MIN_X = 835000;
const double MIN_Y = 815500;
const double INDEX_RANGE_X = -1975;
const double INDEX_RANGE_Y = 43;
// MINH已在query3dRtree.h中定义为宏 -100000

/**
 * @brief 将米转换为度（在给定纬度）
 *
 * 原理：
 * - 纬度1度 ≈ 111,320米（恒定）
 * - 经度1度 ≈ 111,320 * cos(纬度) 米（随纬度变化）
 */
double DEMGenerator::MetersToDegreesAtLatitude(double meters, double latitude) {
    // 纬度方向：1度约等于111,320米
    const double METERS_PER_DEGREE_LAT = 111320.0;

    // 经度方向：随纬度变化，在赤道最大，在极点为0
    // 使用纬度中点的余弦值来计算平均经度分辨率
    double lat_rad = latitude * M_PI / 180.0;
    double meters_per_degree_lon = METERS_PER_DEGREE_LAT * std::cos(lat_rad);

    // 使用经纬度的平均值（简化处理，假设研究区域较小）
    double meters_per_degree_avg = (METERS_PER_DEGREE_LAT + meters_per_degree_lon) / 2.0;

    return meters / meters_per_degree_avg;
}

/**
 * @brief 经纬度（EPSG:4326）转换为R-Tree查询坐标（EPSG:2326）
 *
 * 参考 radar.hpp line 52-76 的实现
 */
bool DEMGenerator::LonLatToRTreeCoord(double lon, double lat, double& rtree_x, double& rtree_y) {
    // 1. 执行PROJ转换（EPSG:4326 → EPSG:2326）
    PJ_CONTEXT* context = proj_context_create();
    PJ* proj_from = proj_create_crs_to_crs(context, "EPSG:4326", "EPSG:2326", nullptr);

    if (!proj_from) {
        std::cerr << "Failed to create PROJ transformer (EPSG:4326 -> EPSG:2326)" << std::endl;
        proj_context_destroy(context);
        return false;
    }

    // 2. 注意：PROJ使用 (lat, lon) 顺序！
    PJ_COORD coord_start = proj_coord(lat, lon, 0, 0);
    PJ_COORD result_start = proj_trans(proj_from, PJ_FWD, coord_start);

    proj_destroy(proj_from);
    proj_context_destroy(context);

    // 3. 计算R-Tree查询坐标（参考 radar.hpp line 71-75）
    // result_start.xy.y 是 X 坐标（东西向）
    // result_start.xy.x 是 Y 坐标（南北向）
    double local_x = result_start.xy.y - MIN_X;
    double local_y = result_start.xy.x - MIN_Y;

    rtree_x = INDEX_RANGE_X + local_x;
    rtree_y = INDEX_RANGE_Y + local_y;

    return true;
}

/**
 * @brief 查询指定经纬度点的高程（使用R-Tree）
 */
float DEMGenerator::QueryHeightAtLonLat(RTree3d* rtree, double lon, double lat) {
    if (!rtree) {
        return -9999.0f;
    }

    // 1. 经纬度转换为R-Tree查询坐标
    double rtree_x, rtree_y;
    if (!LonLatToRTreeCoord(lon, lat, rtree_x, rtree_y)) {
        return -9999.0f;
    }

    // 2. 使用R-Tree查询高度（参考 radar.hpp line 75）
    double height = rtree->Getheight3d(rtree_x, rtree_y, HeightCallback);

    // 3. 检查是否为有效高度
    if (height < MINH) {
        return -9999.0f;
    }

    return static_cast<float>(height);
}

/**
 * @brief 设置GeoTIFF的地理变换参数（EPSG:4326坐标系）
 */
void DEMGenerator::SetGeoTransformWGS84(
    GDALDataset* dataset,
    double min_lon,
    double max_lat,
    double resolution_deg
) {
    // 地理变换参数（度为单位）
    double geotransform[6] = {
        min_lon,         // [0] 左上角经度
        resolution_deg,  // [1] 经度分辨率（度/像素）
        0,               // [2] 旋转参数
        max_lat,         // [3] 左上角纬度
        0,               // [4] 旋转参数
        -resolution_deg  // [5] 纬度分辨率（负值，因为纬度向下递减）
    };
    dataset->SetGeoTransform(geotransform);

    // 设置EPSG:4326投影
    OGRSpatialReference srs;
    srs.importFromEPSG(4326);  // WGS84
    char* wkt = nullptr;
    srs.exportToWkt(&wkt);
    dataset->SetProjection(wkt);
    CPLFree(wkt);
}

/**
 * @brief 从三角面片R-Tree生成DEM（EPSG:4326坐标系）
 *
 * @param resolution_m 分辨率（米），例如1.0表示1米分辨率
 */
bool DEMGenerator::GenerateDEMFromTriangles(
    RTree3d* rtree,
    double lon_range[2],
    double lat_range[2],
    double resolution_m,
    const char* output_path
) {
    if (!rtree) {
        std::cerr << "Error: RTree is null" << std::endl;
        return false;
    }

    // 计算研究区域中心纬度，用于米到度的转换
    double center_lat = (lat_range[0] + lat_range[1]) / 2.0;

    // 将分辨率从米转换为度
    double resolution_deg = MetersToDegreesAtLatitude(resolution_m, center_lat);

    std::cout << "\n========================================" << std::endl;
    std::cout << "Generating DEM from triangulated surface" << std::endl;
    std::cout << "Output CRS: EPSG:4326 (WGS84)" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Output: " << output_path << std::endl;
    std::cout << "Resolution: " << resolution_m << " m/pixel";
    std::cout << " (~" << resolution_deg << " deg/pixel)" << std::endl;
    std::cout << "Lon range: [" << lon_range[0] << ", " << lon_range[1] << "]" << std::endl;
    std::cout << "Lat range: [" << lat_range[0] << ", " << lat_range[1] << "]" << std::endl;

    // 1. 计算栅格尺寸
    int width = static_cast<int>(std::ceil((lon_range[1] - lon_range[0]) / resolution_deg));
    int height = static_cast<int>(std::ceil((lat_range[1] - lat_range[0]) / resolution_deg));

    std::cout << "Raster size: " << width << " x " << height << " pixels" << std::endl;

    // 2. 初始化GDAL
    GDALAllRegister();

    // 3. 创建GeoTIFF数据集
    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GTiff");
    if (!driver) {
        std::cerr << "Error: GTiff driver not available" << std::endl;
        return false;
    }

    GDALDataset* dataset = driver->Create(
        output_path, width, height, 1, GDT_Float32, nullptr
    );
    if (!dataset) {
        std::cerr << "Error: Failed to create output dataset" << std::endl;
        return false;
    }

    // 4. 设置地理变换和投影信息（EPSG:4326）
    SetGeoTransformWGS84(dataset, lon_range[0], lat_range[1], resolution_deg);

    // 5. 获取栅格波段
    GDALRasterBand* band = dataset->GetRasterBand(1);
    band->SetNoDataValue(-9999.0f);

    // 6. 填充高程数据
    std::cout << "Querying heights from R-Tree (EPSG:2326)..." << std::endl;

    float* scanline = new float[width];
    int progress_interval = std::max(1, height / 20);  // 每5%显示一次进度

    for (int row = 0; row < height; row++) {
        // 显示进度
        if (row % progress_interval == 0) {
            std::cout << "Progress: " << (row * 100 / height) << "%" << std::endl;
        }

        for (int col = 0; col < width; col++) {
            // 计算像素中心的经纬度坐标（EPSG:4326）
            double lon = lon_range[0] + col * resolution_deg + resolution_deg / 2.0;
            double lat = lat_range[1] - row * resolution_deg - resolution_deg / 2.0;

            // 查询该经纬度点的地形高度
            // QueryHeightAtLonLat 内部会转换到 EPSG:2326 查询R-Tree
            scanline[col] = QueryHeightAtLonLat(rtree, lon, lat);
        }

        // 写入一行数据到GeoTIFF
        CPLErr err = band->RasterIO(GF_Write, 0, row, width, 1,
                                     scanline, width, 1, GDT_Float32, 0, 0);
        if (err != CE_None) {
            std::cerr << "Error: Failed to write raster data at row " << row << std::endl;
        }
    }

    delete[] scanline;

    // 7. 关闭数据集
    GDALClose(dataset);

    std::cout << "Progress: 100%" << std::endl;
    std::cout << "DEM generation completed: " << output_path << std::endl;
    std::cout << "========================================\n" << std::endl;

    return true;
}

/**
 * @brief 生成多分辨率DEM集合
 *
 * @param resolutions_m 分辨率数组（米），例如 {0.5, 1.0, 5.0, 10.0}
 */
void DEMGenerator::GenerateMultiResolutionDEMs(
    RTree3d* rtree,
    double lon_range[2],
    double lat_range[2],
    const double* resolutions_m,
    int num_resolutions,
    const char* output_prefix
) {
    std::cout << "\n****************************************" << std::endl;
    std::cout << "Generating Multi-Resolution DEM Set" << std::endl;
    std::cout << "****************************************" << std::endl;
    std::cout << "Number of resolutions: " << num_resolutions << std::endl;

    for (int i = 0; i < num_resolutions; i++) {
        double res_m = resolutions_m[i];

        // 生成文件名
        std::string filename = std::string(output_prefix);
        if (res_m < 1.0) {
            // 小于1米：使用小数表示，如 dem_0.5m.tif
            filename += std::to_string(res_m).substr(0, 3) + "m.tif";
        } else {
            // 大于等于1米：使用整数表示，如 dem_1m.tif
            filename += std::to_string(static_cast<int>(res_m)) + "m.tif";
        }

        std::cout << "\n[" << (i+1) << "/" << num_resolutions << "] Generating "
                  << filename << " (resolution: " << res_m << " m)" << std::endl;

        bool success = GenerateDEMFromTriangles(
            rtree, lon_range, lat_range, res_m, filename.c_str()
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
