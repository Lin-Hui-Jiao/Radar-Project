#pragma once

#include <gdal_priv.h>
#include <string>

namespace hiradar {

/**
 * @brief DEM数据加载与查询类
 *
 * 【重要】坐标系统说明：
 * - 输入坐标系必须与DEM文件坐标系一致
 * - 如果DEM是EPSG:4326，则输入 (x=经度, y=纬度) 十进制度
 * - 如果DEM是EPSG:2326，则输入 (x=东向, y=北向) 米
 * - 高程：绝对高程（米，ASL - Above Sea Level）
 *
 * 【功能】：
 * - 加载GeoTIFF格式的DEM文件
 * - 提供高程查询接口（支持双线性插值）
 * - 管理DEM元数据（分辨率、范围等）
 * - 整个DEM加载到内存以提高查询速度
 *
 * 【使用场景】：
 * - 被DEMVisibility使用进行2.5D可见性分析
 * - 可独立使用查询任意点的高程
 */
class DEMLoader {
public:
    DEMLoader();
    ~DEMLoader();

    // 禁止拷贝
    DEMLoader(const DEMLoader&) = delete;
    DEMLoader& operator=(const DEMLoader&) = delete;

    /**
     * @brief 加载GeoTIFF格式的DEM文件
     *
     * @param filepath DEM文件路径
     * @return 是否加载成功
     */
    bool Load(const std::string& filepath);

    /**
     * @brief 查询指定坐标的高程值（双线性插值）
     *
     * 【重要】坐标系要求：
     * - 输入坐标 (x, y) 必须与DEM文件的坐标系一致！
     * - 如果DEM是EPSG:4326: x=经度(十进制度), y=纬度(十进制度)
     * - 如果DEM是EPSG:2326: x=东向距离(米), y=北向距离(米)
     * - 调用者需要自行确保坐标系匹配，本类不做坐标转换
     *
     * 【实现细节】：
     * - 使用像素中心坐标系进行插值（符合GIS标准）
     * - 自动修正GDAL GeoTransform的像素角点偏移
     * - 在边界附近使用最近邻插值兜底
     * - 对NoData像素使用部分有效值平均
     *
     * 【典型用法】：
     * ```cpp
     * // 对于EPSG:4326 DEM
     * double lon = 114.1694;  // 经度（十进制度）
     * double lat = 22.3193;   // 纬度（十进制度）
     * float height = dem.GetHeight(lon, lat);
     * ```
     *
     * @param x X坐标（经度或东向距离，取决于DEM坐标系）
     * @param y Y坐标（纬度或北向距离，取决于DEM坐标系）
     * @return 绝对高程（米，ASL），如果超出范围或全NoData返回nodata_value_
     */
    float GetHeight(double x, double y) const;

    /**
     * @brief 查询指定像素位置的高程（直接读取，无插值）
     *
     * @param col, row 像素坐标
     * @return 高程值（米）
     */
    float GetHeightAtPixel(int col, int row) const;

    /**
     * @brief 检查坐标是否在DEM范围内
     */
    bool IsInBounds(double x, double y) const;

    // 获取DEM元数据
    double GetResolution() const { return resolution_; }
    void GetBounds(double bounds[4]) const;
    int GetWidth() const { return width_; }
    int GetHeight() const { return height_; }
    float GetNoDataValue() const { return nodata_value_; }
    bool IsLoaded() const { return dataset_ != nullptr; }

private:
    GDALDataset* dataset_;
    GDALRasterBand* band_;
    float* data_;           // 整个DEM数据缓存到内存（提高查询速度）
    int width_, height_;
    double geotransform_[6];
    double resolution_;
    float nodata_value_;

    /**
     * @brief 世界坐标 → 像素中心坐标
     *
     * 关键修正：
     * 输出的col和row是相对于像素中心的坐标，而不是像素角点
     * 这确保了双线性插值的正确性，避免了0.5像素的系统性偏移
     *
     * 例如：col=0.0表示查询点正好在像素0的中心
     *       col=0.5表示查询点在像素0和像素1两个中心的正中间
     */
    void WorldToPixel(double x, double y, double& col, double& row) const;

    /**
     * @brief 双线性插值（基于像素中心）
     *
     * 输入：col和row为相对于像素中心的坐标（来自WorldToPixel）
     * 输出：在四个像素中心之间进行双线性插值得到的高程值
     */
    float BilinearInterpolate(double col, double row) const;
};

} // namespace hiradar
