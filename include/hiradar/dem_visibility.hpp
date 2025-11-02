#pragma once

#include "hiradar/dem_loader.hpp"
#include "hiradar/occlusion_utils.h"
#include <vector>

namespace hiradar {

/**
 * @brief 基于DEM的2.5D可见性分析器
 *
 * 【重要】坐标系统说明：
 * - 输入坐标必须为 EPSG:4326 (WGS84 经纬度，十进制度格式)
 * - Vec3d.x: 经度（十进制度）
 * - Vec3d.y: 纬度（十进制度）
 * - Vec3d.z: 绝对高程（米，ASL - Above Sea Level）
 *   注意：必须是绝对高程，不是相对于地面的高度！
 *
 * 【算法说明】2.5D射线采样遮挡检测：
 * 1. 在雷达和目标之间的经纬度空间绘制射线
 * 2. 沿射线均匀采样N个点
 * 3. 对每个采样点查询DEM高度（双线性插值）
 * 4. 如果射线高度 < 地形高度-0.05m，判定为遮挡
 *
 * 【与R-Tree 3D方法的差异】：
 * - R-Tree方法：精确3D几何相交检测（射线-三角面片）
 * - DEM方法：2.5D近似算法（射线采样+DEM查询）
 * - 在复杂地形（山脊、悬崖）下可能产生不同结果
 * - DEM方法更快但精度略低，适合快速可视化和初步分析
 * - 实验对比时差异越大越能体现三角面片方法的优势
 */
class DEMVisibility {
public:
    explicit DEMVisibility(const DEMLoader* dem) : dem_(dem) {}

    /**
     * @brief 判断目标点是否被地形遮挡（2.5D算法）
     *
     * 【重要】输入要求：
     * - radar_pos.x/y: 雷达位置经纬度（十进制度）
     * - radar_pos.z: 雷达绝对高程（米，ASL）
     * - target_pos.x/y: 目标位置经纬度（十进制度）
     * - target_pos.z: 目标绝对高程（米，ASL）
     *
     * 【算法细节】：
     * - 采样密度自动调整：至少1个采样点/DEM像素
     * - 增加2倍安全系数以提高检测精度
     * - 跳过起点和终点（避免自遮挡）
     * - NoData区域跳过（保守策略）
     * - 5cm容差避免浮点精度误判
     *
     * @param radar_pos 雷达位置（Vec3d结构，绝对高程）
     * @param target_pos 目标位置（Vec3d结构，绝对高程）
     * @param sample_count 基础采样点数（默认100，实际会动态调整）
     * @return true=被遮挡（不可见）, false=可见
     */
    bool IsOccluded(const Vec3d& radar_pos,
                    const Vec3d& target_pos,
                    int sample_count = 100) const;

    /**
     * @brief 遮挡分析详细结果（用于调试和可视化）
     */
    struct OcclusionDetail {
        std::vector<Vec3d> sample_points;      // 采样点位置（十进制度+绝对高程）
        std::vector<float> terrain_heights;    // 地形高度（米，从DEM查询）
        std::vector<float> ray_heights;        // 射线高度（米，绝对高程）
        std::vector<bool> is_below_terrain;    // 是否在地形下方
        bool final_occluded;                   // 最终判断结果
        int first_occlusion_index;             // 第一个遮挡点索引（-1表示无遮挡）
        double horizontal_distance_m;          // 水平距离（米）
        int actual_sample_count;               // 实际采样点数
    };

    /**
     * @brief 执行遮挡分析并返回详细信息（用于调试）
     *
     * 返回完整的采样点信息，包括每个点的：
     * - 位置坐标
     * - 地形高度
     * - 射线高度
     * - 是否被遮挡
     *
     * @param radar_pos 雷达位置
     * @param target_pos 目标位置
     * @param sample_count 基础采样点数
     * @return 详细的遮挡分析结果
     */
    OcclusionDetail AnalyzeOcclusion(const Vec3d& radar_pos,
                                     const Vec3d& target_pos,
                                     int sample_count = 100) const;

private:
    const DEMLoader* dem_;

    /**
     * @brief 沿射线均匀采样N个点（内部辅助函数）
     *
     * 在经纬度空间进行线性插值（近似测地线）
     * 注意：在小范围（<50km）误差可忽略
     *
     * @param start 起点（十进制度+绝对高程）
     * @param end 终点（十进制度+绝对高程）
     * @param count 采样点数（不含起点终点）
     * @return 采样点列表（包含起点和终点，共count+1个点）
     */
    std::vector<Vec3d> SampleRayPoints(const Vec3d& start,
                                       const Vec3d& end,
                                       int count) const;
};

} // namespace hiradar
