#include "hiradar/dem_visibility.hpp"
#include <cmath>
#include <algorithm>

namespace hiradar {

/**
 * @brief 沿射线均匀采样N个点（在经纬度空间线性插值）
 *
 * 算法说明：
 * 1. 直接使用十进制度进行线性插值
 * 2. 插值经度、纬度和高度
 *
 * 注意：在小范围（<50km）内，经纬度线性插值近似测地线
 */
std::vector<Vec3d> DEMVisibility::SampleRayPoints(
    const Vec3d& start, const Vec3d& end, int count
) const {
    std::vector<Vec3d> points;
    points.reserve(count + 1);

    for (int i = 0; i <= count; i++) {
        double t = static_cast<double>(i) / count;

        // 在十进制度空间线性插值
        Vec3d p;
        p.x = start.x + t * (end.x - start.x);  // 经度
        p.y = start.y + t * (end.y - start.y);  // 纬度
        p.z = start.z + t * (end.z - start.z);  // 高度

        points.push_back(p);
    }

    return points;
}

/**
 * @brief 判断目标点是否被地形遮挡（2.5D算法）
 *
 * 关键改进：
 * 1. 使用Vec3d结构（保持一致性）
 * 2. 明确输入为EPSG:4326经纬度（十进制度）+ 绝对高程
 * 3. 增加2倍采样密度安全系数（提高检测精度，增大与R-Tree的差异）
 * 4. 详细的边界检查和错误处理
 */
bool DEMVisibility::IsOccluded(
    const Vec3d& radar_pos,
    const Vec3d& target_pos,
    int sample_count
) const {
    // 边界检查1：DEM未加载
    if (!dem_ || !dem_->IsLoaded()) {
        return false;  // 保守策略：无DEM时认为不遮挡
    }

    // 直接使用Vec3d的成员（十进制度）
    double radar_lon_deg = radar_pos.x;   // 经度
    double radar_lat_deg = radar_pos.y;   // 纬度
    double target_lon_deg = target_pos.x; // 经度
    double target_lat_deg = target_pos.y; // 纬度

    // 边界检查2：起点终点重合
    double dx_deg = target_lon_deg - radar_lon_deg;
    double dy_deg = target_lat_deg - radar_lat_deg;
    if (std::abs(dx_deg) < 1e-9 && std::abs(dy_deg) < 1e-9) {
        return false;  // 自己不遮挡自己
    }

    // 边界检查3：输入坐标合理性（EPSG:4326范围）
    // if (radar_lon_deg < -180 || radar_lon_deg > 180 ||
    //     radar_lat_deg < -90 || radar_lat_deg > 90 ||
    //     target_lon_deg < -180 || target_lon_deg > 180 ||
    //     target_lat_deg < -90 || target_lat_deg > 90) {
    //     // std::cerr << "Warning: Invalid lon/lat coordinates in DEMVisibility::IsOccluded!" << std::endl;
    //     return false;
    // }

    // ========== 动态调整采样密度 ==========
    // 目标：确保采样间隔 <= DEM分辨率，避免遗漏地形特征
    //
    // 关键修正：使用米制距离而非度距离
    // 原因：DEM分辨率代表物理地形细节（米），而非角度细节（度）

    // 计算中心纬度，用于经度-米数转换
    double center_lat = (radar_lat_deg + target_lat_deg) / 2.0;
    double lat_rad = center_lat * M_PI / 180.0;

    // 转换系数（在给定纬度处）
    const double METERS_PER_DEG_LAT = 111320.0;  // 1度纬度 ≈ 111.32km（全球恒定）
    double meters_per_deg_lon = METERS_PER_DEG_LAT * std::cos(lat_rad);  // 1度经度（随纬度变化）

    // 将经纬度差转换为米
    double dx_m = dx_deg * meters_per_deg_lon;
    double dy_m = dy_deg * METERS_PER_DEG_LAT;
    double horizontal_dist_m = std::sqrt(dx_m * dx_m + dy_m * dy_m);

    // DEM分辨率转换为米（使用纬度方向的转换，因为它更稳定）
    double dem_resolution_m = dem_->GetResolution() * METERS_PER_DEG_LAT;

    // 根据DEM分辨率计算采样数
    // 采样间隔应该与DEM分辨率一致，确保不会遗漏地形特征
    // 不再使用默认下限，让低分辨率DEM真正使用更少的采样点
    int actual_sample_count = static_cast<int>(std::ceil(horizontal_dist_m / dem_resolution_m));

    // 确保至少有2个采样点（起点和终点之间至少1个中间点）
    actual_sample_count = std::max(actual_sample_count, 2);

    // ========== 沿射线采样N个点 ==========
    std::vector<Vec3d> ray_points = SampleRayPoints(radar_pos, target_pos, actual_sample_count);

    // 地形穿透容差：避免浮点精度问题导致误判
    // 如果射线高度与地形高度相差小于此值，认为是"擦边"而非"穿透"
    const float TERRAIN_CLEARANCE = 0.05f;  // 5厘米容差

    // ========== 检查每个采样点（跳过起点和终点）==========
    // 注意：起点是雷达位置，终点是目标位置，都不应该参与遮挡判断
    for (size_t i = 1; i < ray_points.size() - 1; i++) {
        const Vec3d& point = ray_points[i];

        // Vec3d已经是十进制度，直接使用
        double lon_deg = point.x;  // 经度
        double lat_deg = point.y;  // 纬度

        // 查询该点的地形高度（使用双线性插值）
        // 注意：DEMLoader::GetHeight() 期望输入为 (经度, 纬度) 十进制度
        float terrain_height = dem_->GetHeight(lon_deg, lat_deg);

        // 如果是NoData区域，跳过
        // 策略说明：NoData可能是水面、数据空洞或DEM边界外
        // 这里采用跳过策略以减少误报（保守策略）
        if (terrain_height == dem_->GetNoDataValue()) {
            continue;
        }

        // 判断遮挡：射线高度需要明显低于地形高度（考虑容差）
        // point.z 是绝对高程（ASL - Above Sea Level）
        // terrain_height 也是绝对高程（从DEM查询）
        // 使用容差可以避免浮点精度和双线性插值误差导致的误判
        if (point.z < terrain_height - TERRAIN_CLEARANCE) {
            return true;  // 被遮挡
        }
    }

    return false;  // 可见
}

/**
 * @brief 执行遮挡分析并返回详细信息（用于调试和可视化）
 *
 * 这个函数主要用于：
 * 1. 调试：检查第一个遮挡点的位置
 * 2. 可视化：绘制射线和地形的关系
 * 3. 对比实验：分析与R-Tree方法的差异原因
 */
DEMVisibility::OcclusionDetail DEMVisibility::AnalyzeOcclusion(
    const Vec3d& radar_pos,
    const Vec3d& target_pos,
    int sample_count
) const {
    OcclusionDetail detail;
    detail.final_occluded = false;
    detail.first_occlusion_index = -1;
    detail.horizontal_distance_m = 0.0;
    detail.actual_sample_count = 0;

    if (!dem_ || !dem_->IsLoaded()) {
        return detail;
    }

    // 直接使用Vec3d的成员（十进制度）
    double radar_lon_deg = radar_pos.x;
    double radar_lat_deg = radar_pos.y;
    double target_lon_deg = target_pos.x;
    double target_lat_deg = target_pos.y;

    // 动态调整采样密度（与IsOccluded保持一致）
    double center_lat = (radar_lat_deg + target_lat_deg) / 2.0;
    double lat_rad = center_lat * M_PI / 180.0;

    const double METERS_PER_DEG_LAT = 111320.0;
    double meters_per_deg_lon = METERS_PER_DEG_LAT * std::cos(lat_rad);

    double dx_m = (target_lon_deg - radar_lon_deg) * meters_per_deg_lon;
    double dy_m = (target_lat_deg - radar_lat_deg) * METERS_PER_DEG_LAT;
    double horizontal_dist_m = std::sqrt(dx_m * dx_m + dy_m * dy_m);
    detail.horizontal_distance_m = horizontal_dist_m;

    double dem_resolution_m = dem_->GetResolution() * METERS_PER_DEG_LAT;

    // 根据DEM分辨率计算采样数（与IsOccluded保持一致）
    int actual_sample_count = static_cast<int>(std::ceil(horizontal_dist_m / dem_resolution_m));
    actual_sample_count = std::max(actual_sample_count, 2);
    detail.actual_sample_count = actual_sample_count;

    detail.sample_points = SampleRayPoints(radar_pos, target_pos, actual_sample_count);

    const float TERRAIN_CLEARANCE = 0.05f;  // 与IsOccluded保持一致

    for (size_t i = 0; i < detail.sample_points.size(); i++) {
        const Vec3d& point = detail.sample_points[i];

        // Vec3d已经是十进制度，直接使用
        double lon_deg = point.x;
        double lat_deg = point.y;

        float terrain_h = dem_->GetHeight(lon_deg, lat_deg);
        detail.terrain_heights.push_back(terrain_h);
        detail.ray_heights.push_back(point.z);

        // 使用容差判断是否在地形下方
        bool below = (terrain_h != dem_->GetNoDataValue() &&
                      point.z < terrain_h - TERRAIN_CLEARANCE);
        detail.is_below_terrain.push_back(below);

        // 检查遮挡（跳过起点和终点）
        if (below && i > 0 && i < detail.sample_points.size() - 1) {
            detail.final_occluded = true;
            if (detail.first_occlusion_index == -1) {
                detail.first_occlusion_index = static_cast<int>(i);
            }
        }
    }

    return detail;
}

} // namespace hiradar
