#pragma once

#include <gdal_priv.h>
#include <ogr_spatialref.h>
#include <string>
#include "hiradar/RTree.h"
#include "hiradar/occlusion_utils.h"

namespace hiradar {

/**
 * @brief DEM生成器 - 从三角面片数据栅格化生成DEM
 *
 * 【坐标系统说明】：
 * - 输出：EPSG:4326 (WGS84 经纬度坐标系)
 * - 查询：自动转换到EPSG:2326（香港投影坐标系）查询R-Tree
 * - 高程：绝对高程（米，ASL - Above Sea Level）
 *
 * 【生成流程】：
 * 1. 指定经纬度范围（EPSG:4326，十进制度）和分辨率（米）
 * 2. 生成规则网格（按DEM像素分辨率）
 * 3. 对每个网格点：
 *    a. 经纬度(EPSG:4326) → PROJ转换 → 投影坐标(EPSG:2326)
 *    b. 应用局部偏移 → R-Tree查询坐标
 *    c. 查询三角面片高度 → 绝对高程
 * 4. 输出为GeoTIFF格式DEM（EPSG:4326坐标系）
 *
 * 【用途】：
 * - 用于2.5D可见性分析（DEMVisibility）
 * - 用于可视化和快速预览
 * - 用于与精确3D方法（R-Tree）对比实验
 */
class DEMGenerator {
public:
    /**
     * @brief 从三角面片R-Tree生成DEM（输出EPSG:4326坐标系）
     *
     * 【重要】坐标系统说明：
     * - 输入：经纬度范围（EPSG:4326，十进制度）
     * - 输出：GeoTIFF文件，EPSG:4326坐标系
     * - 查询：内部自动转换到EPSG:2326查询R-Tree三角面片
     * - 高程：绝对高程（米，ASL）
     *
     * 【分辨率说明】：
     * - resolution_m: 物理分辨率（米/像素）
     * - 内部自动转换为角度分辨率（度/像素）
     * - 转换基于研究区域中心纬度
     *
     * 【生成的DEM用途】：
     * - 用于DEMVisibility进行2.5D可见性分析
     * - 用于DEMLoader加载和查询高程
     * - 用于与R-Tree 3D方法对比实验
     *
     * @param rtree 三角面片R-Tree索引（EPSG:2326坐标系）
     * @param lon_range 经度范围 [min_lon, max_lon]（十进制度，EPSG:4326）
     * @param lat_range 纬度范围 [min_lat, max_lat]（十进制度，EPSG:4326）
     * @param resolution_m DEM分辨率（米/像素），例如1.0表示1米分辨率
     * @param output_path 输出GeoTIFF文件路径
     * @return 是否成功生成
     */
    static bool GenerateDEMFromTriangles(
        RTree3d* rtree,
        double lon_range[2],
        double lat_range[2],
        double resolution_m,
        const char* output_path
    );

    /**
     * @brief 生成多分辨率DEM集合
     *
     * @param rtree 三角面片R-Tree索引
     * @param lon_range 经度范围
     * @param lat_range 纬度范围
     * @param resolutions_m 分辨率数组（米），例如 {0.5, 1.0, 5.0, 10.0}
     * @param num_resolutions 分辨率数量
     * @param output_prefix 输出文件名前缀（如"dem_"）
     */
    static void GenerateMultiResolutionDEMs(
        RTree3d* rtree,
        double lon_range[2],
        double lat_range[2],
        const double* resolutions_m,
        int num_resolutions,
        const char* output_prefix
    );

private:
    /**
     * @brief 经纬度（EPSG:4326）转换为R-Tree查询坐标（EPSG:2326）
     *
     * 参考 radar.hpp line 52-76 的实现
     */
    static bool LonLatToRTreeCoord(
        double lon, double lat,
        double& rtree_x, double& rtree_y
    );

    /**
     * @brief 设置GeoTIFF的地理变换参数和投影信息（EPSG:4326）
     *
     * @param resolution_deg 分辨率（度/像素）
     */
    static void SetGeoTransformWGS84(
        GDALDataset* dataset,
        double min_lon, double max_lat,
        double resolution_deg
    );

    /**
     * @brief 查询指定经纬度点的高程（使用R-Tree）
     *
     * 内部会将经纬度转换为EPSG:2326坐标查询R-Tree
     */
    static float QueryHeightAtLonLat(
        RTree3d* rtree,
        double lon, double lat
    );

    /**
     * @brief 将米转换为度（在给定纬度）
     *
     * @param meters 距离（米）
     * @param latitude 纬度（度）
     * @return 对应的度数
     */
    static double MetersToDegreesAtLatitude(double meters, double latitude);
};

} // namespace hiradar
