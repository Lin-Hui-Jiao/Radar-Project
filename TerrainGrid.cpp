#include "hiradar/TerrainGrid.hpp"
#include <iostream>
#include <cmath>
#include <gdal_priv.h> // 引入GDAL头文件
const double PI = 3.14159265358979323846;
const double EARTH_RADIUS_M = 6371000.0;
bool TerrainGrid::build(const std::string& geotiff_path) {
   
    GDALDataset* dataset = (GDALDataset*)GDALOpen(geotiff_path.c_str(), GA_ReadOnly);
    if (!dataset) {
        std::cerr << "Error: Could not open GeoTIFF file " << geotiff_path << std::endl;
        return false;
    }

    // --- 第一步：读取元数据 ---
    this->grid_size_x = dataset->GetRasterXSize();
    this->grid_size_y = dataset->GetRasterYSize();
    
    double geo_transform[6];
    if (dataset->GetGeoTransform(geo_transform) != CE_None) {
        std::cerr << "Error: Could not get geotransform from file." << std::endl;
        GDALClose(dataset);
        return false;
    }
    
    // --- 第二步：定义网格参数 ---
    this->min_lon = geo_transform[0];
    this->cell_width = geo_transform[1];
    this->max_lat = geo_transform[3];
    this->cell_height = geo_transform[5]; // 通常是负值
    //this->min_lat = this->max_lat + this->grid_size_y * this->cell_height;

    // --- 第三步：填充网格 ---
    this->grid.resize(grid_size_x, std::vector<TerrainPoint>(grid_size_y));
    GDALRasterBand* band = dataset->GetRasterBand(1);
    
    float* row_buffer = new float[grid_size_x];

    for (int j = 0; j < grid_size_y; ++j) { // 遍历每一行 (line)
        band->RasterIO(GF_Read, 0, j, grid_size_x, 1, row_buffer, grid_size_x, 1, GDT_Float32, 0, 0); //读取一整行数据到row_buffer中
        for (int i = 0; i < grid_size_x; ++i) { // 遍历每一列 (pixel)
            TerrainPoint point;
            point.alt = row_buffer[i];
            // 计算每个像素中心的经纬度
            point.lon = this->min_lon + (i + 0.5) * this->cell_width;
            point.lat = this->max_lat + (j + 0.5) * this->cell_height;
            
            // 检查计算出的索引是否在范围内
            if (i >= 0 && i < grid_size_x && j >= 0 && j < grid_size_y) {
                 this->grid[i][j] = point;
            }
        }
    }

    delete[] row_buffer;
    GDALClose(dataset);
    std::cout << "Terrain grid built successfully (" << grid_size_x << "x" << grid_size_y << ")." << std::endl;
    return true;
}

float TerrainGrid::getNearestElevation(float lon, float lat) const {
    if (grid.empty()) return -9999.0f; // 未初始化
    // --- 步骤A：定位中心格子 ---
    int i = static_cast<int>((lon - min_lon) / cell_width);
    int j = static_cast<int>((lat - max_lat) / cell_height); // 正确的纬度索引计算

    TerrainPoint best_point = {};
    float min_dist_sq = std::numeric_limits<float>::max();

    // --- 步骤B & C：自适应搜索3x3邻域 ---
    for (int search_j = j - 1; search_j <= j + 1; ++search_j) {
        for (int search_i = i - 1; search_i <= i + 1; ++search_i) {
            // 边界检查
            if (search_i < 0 || search_i >= grid_size_x || search_j < 0 || search_j >= grid_size_y) {
                continue;
            }
            const auto& point = grid[search_i][search_j];
            float dx = lon - point.lon;
            float dy = lat - point.lat;
            float dist_sq = dx * dx + dy * dy;
            if (dist_sq < min_dist_sq) {
                min_dist_sq = dist_sq;
                best_point = point;
            }
        }
    }

    return best_point.alt;
}

// 实现 getGroundResolution 函数
float TerrainGrid::getGroundResolution() const {
    if (grid.empty()) return 1000.0f; // 如果未初始化，返回一个默认值

    // 纬度转弧度
    float min_lat = this->max_lat + this->grid_size_y * this->cell_height;
    double mid_lat_rad = (min_lat + (grid_size_y * cell_height / 2.0)) * PI / 180.0;
    
    // 计算一度经线在当前纬度的长度（米）
    double lon_dist_per_deg = (PI / 180.0) * EARTH_RADIUS_M * cos(mid_lat_rad);
    // 计算一度纬线的长度（米）
    double lat_dist_per_deg = (PI / 180.0) * EARTH_RADIUS_M;

    // 将单元格的度数尺寸转换为米
    float resolution_x_m = std::abs(cell_width * lon_dist_per_deg);
    float resolution_y_m = std::abs(cell_height * lat_dist_per_deg);

    // 返回两者中较小（更精细）的一个作为保守的采样间隔
    return std::min(resolution_x_m, resolution_y_m);
}