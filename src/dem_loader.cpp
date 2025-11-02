#include "hiradar/dem_loader.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>

namespace hiradar {

DEMLoader::DEMLoader()
    : dataset_(nullptr)
    , band_(nullptr)
    , data_(nullptr)
    , width_(0)
    , height_(0)
    , resolution_(0.0)
    , nodata_value_(-9999.0f)
{
    std::fill(geotransform_, geotransform_ + 6, 0.0);
}

DEMLoader::~DEMLoader() {
    if (data_) {
        delete[] data_;
        data_ = nullptr;
    }
    if (dataset_) {
        GDALClose(dataset_);
        dataset_ = nullptr;
    }
}

bool DEMLoader::Load(const std::string& filepath) {
    GDALAllRegister();

    dataset_ = (GDALDataset*)GDALOpen(filepath.c_str(), GA_ReadOnly);
    if (!dataset_) {
        std::cerr << "Error: Cannot open DEM file: " << filepath << std::endl;
        return false;
    }

    // 读取元数据
    width_ = dataset_->GetRasterXSize();
    height_ = dataset_->GetRasterYSize();
    dataset_->GetGeoTransform(geotransform_);
    resolution_ = geotransform_[1];  // 假设X和Y分辨率相同

    band_ = dataset_->GetRasterBand(1);
    if (!band_) {
        std::cerr << "Error: Cannot access raster band" << std::endl;
        GDALClose(dataset_);
        dataset_ = nullptr;
        return false;
    }

    int hasNoData = 0;
    nodata_value_ = static_cast<float>(band_->GetNoDataValue(&hasNoData));
    if (!hasNoData) {
        nodata_value_ = -9999.0f;  // 默认NoData值
    }

    // 将整个DEM加载到内存（提高查询速度）
    data_ = new float[width_ * height_];
    CPLErr err = band_->RasterIO(GF_Read, 0, 0, width_, height_,
                                  data_, width_, height_, GDT_Float32, 0, 0);

    if (err != CE_None) {
        std::cerr << "Error: Failed to read DEM data" << std::endl;
        delete[] data_;
        data_ = nullptr;
        GDALClose(dataset_);
        dataset_ = nullptr;
        return false;
    }

    std::cout << "DEM loaded successfully: " << filepath << std::endl;
    std::cout << "  Size: " << width_ << " x " << height_ << " pixels" << std::endl;
    std::cout << "  Resolution: " << resolution_ << " m/pixel" << std::endl;

    double bounds[4];
    GetBounds(bounds);
    std::cout << "  Bounds: [" << bounds[0] << ", " << bounds[1]
              << ", " << bounds[2] << ", " << bounds[3] << "]" << std::endl;

    return true;
}

void DEMLoader::WorldToPixel(double x, double y, double& col, double& row) const {
    // 应用仿射变换的逆运算，转换到像素中心坐标系
    //
    // GDAL GeoTransform约定：
    //   geotransform[0] 和 geotransform[3] 是左上角像素的左上角点坐标
    //
    // 像素中心坐标：
    //   像素[i,j]的中心在世界坐标系中的位置为：
    //   x_center = geotransform[0] + (i + 0.5) * geotransform[1]
    //   y_center = geotransform[3] + (j + 0.5) * geotransform[5]
    //
    // 因此，从世界坐标反推像素索引（相对于像素中心）：
    //   i = (x - geotransform[0]) / geotransform[1] - 0.5
    //   j = (y - geotransform[3]) / geotransform[5] - 0.5
    //
    // 这样，当查询点正好位于像素中心时，col和row为整数
    // 当查询点位于四个像素中心围成的区域内时，col和row为[n, n+1)之间的小数

    col = (x - geotransform_[0]) / geotransform_[1] - 0.5;
    row = (y - geotransform_[3]) / geotransform_[5] - 0.5;
}

float DEMLoader::BilinearInterpolate(double col, double row) const {
    // 双线性插值：在四个像素中心之间插值
    //
    // 输入的col和row已经是相对于像素中心的坐标（来自WorldToPixel）
    // 例如：col=0.3表示查询点位于像素0中心右侧0.3个像素宽度处
    //
    // 我们需要找到包围查询点的四个像素中心：
    //   col0 = floor(col)     // 左侧像素索引
    //   col1 = col0 + 1       // 右侧像素索引
    //   row0 = floor(row)     // 上侧像素索引
    //   row1 = row0 + 1       // 下侧像素索引

    int col0 = static_cast<int>(std::floor(col));
    int row0 = static_cast<int>(std::floor(row));
    int col1 = col0 + 1;
    int row1 = row0 + 1;

    // 边界检查：确保四个像素都在有效范围内
    // 注意：col0和row0可能为负数（查询点在第一个像素中心左侧或上侧）
    if (col0 < 0 || row0 < 0 || col1 >= width_ || row1 >= height_) {
        // 如果任何一个像素超出范围，尝试使用最近邻插值
        // 这样可以在边界附近也能获得合理的高程值
        int nearest_col = std::max(0, std::min(width_ - 1, static_cast<int>(std::round(col))));
        int nearest_row = std::max(0, std::min(height_ - 1, static_cast<int>(std::round(row))));

        float h = data_[nearest_row * width_ + nearest_col];
        return (h == nodata_value_) ? nodata_value_ : h;
    }

    // 获取四个像素中心的高程值
    float h00 = data_[row0 * width_ + col0];  // 左上
    float h10 = data_[row0 * width_ + col1];  // 右上
    float h01 = data_[row1 * width_ + col0];  // 左下
    float h11 = data_[row1 * width_ + col1];  // 右下

    // 检查NoData值
    if (h00 == nodata_value_ || h10 == nodata_value_ ||
        h01 == nodata_value_ || h11 == nodata_value_) {
        // 如果任何一个像素是NoData，尝试使用有效值的平均值
        int valid_count = 0;
        float sum = 0.0f;
        if (h00 != nodata_value_) { sum += h00; valid_count++; }
        if (h10 != nodata_value_) { sum += h10; valid_count++; }
        if (h01 != nodata_value_) { sum += h01; valid_count++; }
        if (h11 != nodata_value_) { sum += h11; valid_count++; }

        if (valid_count == 0) {
            return nodata_value_;
        }
        // 至少有一个有效值时，返回平均值而不是NoData
        return sum / valid_count;
    }

    // 计算插值权重（查询点在四个像素中心围成的矩形内的相对位置）
    // dx, dy ∈ [0, 1]
    double dx = col - col0;  // 相对于左侧像素中心的横向偏移比例
    double dy = row - row0;  // 相对于上侧像素中心的纵向偏移比例

    // 双线性插值公式：
    // h(x,y) = h00*(1-dx)*(1-dy) + h10*dx*(1-dy) + h01*(1-dx)*dy + h11*dx*dy
    //
    // 含义：按照查询点到四个像素中心的距离进行加权平均
    float h = static_cast<float>(
        h00 * (1.0 - dx) * (1.0 - dy) +  // 左上角权重
        h10 * dx * (1.0 - dy) +          // 右上角权重
        h01 * (1.0 - dx) * dy +          // 左下角权重
        h11 * dx * dy                    // 右下角权重
    );

    return h;
}

float DEMLoader::GetHeight(double x, double y) const {
    if (!data_) {
        return nodata_value_;
    }

    double col, row;
    WorldToPixel(x, y, col, row);
    return BilinearInterpolate(col, row);
}

float DEMLoader::GetHeightAtPixel(int col, int row) const {
    if (!data_ || col < 0 || row < 0 || col >= width_ || row >= height_) {
        return nodata_value_;
    }
    return data_[row * width_ + col];
}
//这是一个工具函数，用于快速判断一个世界坐标 (x, y) 是否落在整个 DEM 文件的地理范围之内。
bool DEMLoader::IsInBounds(double x, double y) const {
    double col, row;
    WorldToPixel(x, y, col, row);
    return (col >= 0 && col < width_ && row >= 0 && row < height_);
}
//获取并返回 DEM 数据覆盖的完整地理范围，即边界框（Bounding Box）的四个角点坐标 [minX, minY, maxX, maxY]。
void DEMLoader::GetBounds(double bounds[4]) const {
    // bounds[0] = minX, bounds[1] = minY, bounds[2] = maxX, bounds[3] = maxY
    bounds[0] = geotransform_[0];  // 左上角X
    bounds[1] = geotransform_[3] + height_ * geotransform_[5];  // 左下角Y（geotransform_[5]是负值）
    bounds[2] = geotransform_[0] + width_ * geotransform_[1];   // 右上角X
    bounds[3] = geotransform_[3];  // 左上角Y
}

} // namespace hiradar
