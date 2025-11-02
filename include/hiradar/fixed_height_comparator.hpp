#pragma once

#include "hiradar/RTree.h"
#include "hiradar/dem_loader.hpp"
#include "hiradar/dem_visibility.hpp"
#include "hiradar/occlusion_utils.h"
#include <string>
#include <vector>

// 前向声明，避免循环依赖
class Radar;

namespace hiradar {

/**
 * @brief 固定高度平面的可见性和能量密度对比实验
 *
 * 功能说明：
 * 1. 在固定高度平面上以1m分辨率均匀采样网格点
 * 2. 对每个网格点：
 *    - 使用三角面片判断可见性
 *    - 使用DEM数据判断可见性
 *    - 记录可见性差异（disagreement）
 * 3. 对多个不同分辨率的DEM重复实验
 * 4. 生成PNG可视化结果，使用颜色编码表示：
 *    - 灰色：被遮挡（两种方法一致）
 *    - 彩色：雷达能量密度（三角面片判断可见）
 *    - 特殊标记：可见性差异点（两种方法不一致）
 */
class FixedHeightComparator {
public:
    /**
     * @brief 网格点信息结构
     */
    struct GridPoint {
        Vec3d position;        // 位置（EPSG:4326: x=lon, y=lat, z=alt）
        bool visible_mesh;     // 三角面片判断：是否可见
        bool visible_dem;      // DEM判断：是否可见
        float power_density;   // 雷达能量密度（仅当mesh判断可见时有效）
    };

    /**
     * @brief 实验结果结构
     */
    struct ExperimentResult {
        std::vector<GridPoint> grid_points;  // 所有网格点

        // 统计信息
        size_t total_points;                 // 总点数
        size_t visible_mesh_count;           // 三角面片判断可见数
        size_t visible_dem_count;            // DEM判断可见数
        size_t disagreement_count;           // 差异点数

        // 能量密度范围（用于颜色映射）
        float min_power_density;
        float max_power_density;

        // 实验参数
        double fixed_height;                 // 固定高度（米）
        double grid_resolution;              // 网格分辨率（米）
        std::string dem_file;                // 使用的DEM文件
        double dem_resolution;               // DEM分辨率（度）

        // 性能统计（时间单位：秒）
        double total_time;                   // 总耗时
        double mesh_visibility_time;         // 三角面片可见性判断耗时
        double dem_visibility_time;          // DEM可见性判断耗时
        double power_calculation_time;       // 能量密度计算耗时
        int num_threads;                     // 使用的OpenMP线程数
    };

    /**
     * @brief 执行固定高度平面对比实验
     *
     * @param radar 雷达对象指针（用于调用能量密度计算接口）
     * @param radar_pos 雷达位置（EPSG:4326: x=lon°, y=lat°, z=alt_m，其中z为绝对海拔高度）
     * @param lon_range 经度范围 [min_lon, max_lon]（度）
     * @param lat_range 纬度范围 [min_lat, max_lat]（度）
     * @param fixed_height 固定高度（米，绝对海拔高度ASL - Above Sea Level）
     * @param grid_resolution_m 网格采样分辨率（米，建议1m）
     * @param rtree 三角面片R-Tree索引
     * @param dem DEM数据加载器
     * @return ExperimentResult 实验结果
     */
    static ExperimentResult RunExperiment(
        Radar* radar,
        const Vec3d& radar_pos,
        double lon_range[2],
        double lat_range[2],
        double fixed_height,
        double grid_resolution_m,
        RTree3d* rtree,
        DEMLoader* dem
    );

    /**
     * @brief 生成PNG可视化结果
     *
     * 颜色编码规则：
     * - 灰色 (128,128,128): 被遮挡（两种方法一致判断为不可见）
     * - 黄-红渐变: 雷达能量密度（三角面片判断可见）
     *   - 浅黄色：最低能量
     *   - 纯红色：最高能量
     * - 蓝色 (0,0,255): 仅DEM判断可见，三角面片判断不可见（假阴性）
     * - 绿色 (0,255,0): 仅三角面片判断可见，DEM判断不可见（假阳性）
     *
     * @param result 实验结果
     * @param output_file 输出PNG文件路径
     */
    static void GenerateVisualization(
        const ExperimentResult& result,
        const std::string& output_file
    );

    /**
     * @brief 导出详细CSV结果
     *
     * CSV格式：
     * lon,lat,alt,visible_mesh,visible_dem,disagreement,power_density
     *
     * @param result 实验结果
     * @param output_file 输出CSV文件路径
     */
    static void ExportResultsToCSV(
        const ExperimentResult& result,
        const std::string& output_file
    );

    /**
     * @brief 打印结果摘要
     */
    static void PrintSummary(const ExperimentResult& result);

private:
    /**
     * @brief 将经纬度距离转换为米
     *
     * @param degrees 度数
     * @param latitude 中心纬度（用于经度转换）
     * @param is_longitude true=经度差，false=纬度差
     * @return 米数
     */
    static double DegreesToMeters(double degrees, double latitude, bool is_longitude);

    /**
     * @brief 将米距离转换为经纬度
     *
     * @param meters 米数
     * @param latitude 中心纬度（用于经度转换）
     * @param is_longitude true=经度，false=纬度
     * @return 度数
     */
    static double MetersToDegrees(double meters, double latitude, bool is_longitude);

    // 注意：原有的CalculatePowerDensity函数已被移除
    // 现在统一使用Radar类的CalculateSinglePointPowerDensityFromVec3d接口
};

} // namespace hiradar
