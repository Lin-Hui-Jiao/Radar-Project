#pragma once
#include <map>
#include "hiradar/interface.hpp"
#include<cmath>
#include<vector>
#include "occlusion_utils.h" 
static const float LevelStepMap[][3] = {
    {256, 0, 0},    //0级
    {256, 0, 0},    //1级
    {128, 0, 0},
    {64, 0, 0},
    {32, 0, 0},
    {16, 0, 0},
    {8, 0, 0},
    {4, 0, 0},
    {2, 0, 0},
    {1, 0, 0},      //9级，0到9级的划分都是以度为单位
    {0, 32, 0},     //10级
    {0, 16, 0},
    {0, 8, 0},
    {0, 4, 0},
    {0, 2, 0},
    {0, 1, 0},      //15级，10到15级的划分都是以分为单位
    {0, 0, 32},     //16级
    {0, 0, 16},
    {0, 0, 8},
    {0, 0, 4},
    {0, 0, 2},
    {0, 0, 1}};     //21级，16到21级的划分都是以秒为单位


inline void DecimalToDMS(float dms[3], float de)
{
     // 度：直接取整数部分
     int D = static_cast<int>(de);
    
     // 分：用小数部分乘以60，然后取整数部分
     int M = static_cast<int>((de - D) * 60);
     
     // 秒：用度的小数部分乘以3600，减去已经计算过的分的整数秒数
     // (de - D) * 3600 是总秒数
     // M * 60 是“分”所占的秒数
     float S = (de - D) * 3600.0f - M * 60.0f;
     
     dms[0] = D;
     dms[1] = M;
     dms[2] = S;
}

inline float DMSLessThan(float dms[3], float other[3])
{
    return (dms[0] < other[0]) ||
           (dms[0] == other[0] &&
            (dms[1] < other[1] || (dms[1] == other[1] && dms[2] < other[2])));
}

//高程（单位为米）转化为度分秒形式
inline void AltToDMS(float dms[3],float alt)
{
    float alt_dec=((alt/1000)+EARTH_RADIUS)/110;
    int D=alt_dec;
    int F=(alt_dec-D)*60;
    int E=((alt_dec-D)*60-F)*60;
    dms[0]=D;
    dms[1]=F;
    dms[2]=E;
}

//度秒分转化为高程,单位为米
inline float DMSToAlt(const float dms[3])
{
    float alt_dec = dms[0] + dms[1] / 60 + dms[2] / 3600;
    return (alt_dec*110-EARTH_RADIUS)*1000;
}



//0到9级，以度作步长
inline std::pair<Position*, int> CreateGridbyDe(Position radcenter, float range, unsigned short level)
{
    float radlat_de = DMSToDecimal(radcenter.lat);
    float radlon_de = DMSToDecimal(radcenter.lon);
    //先将radcenter的高程转化为度分秒的形式
    float radalt_dms[3] = { 0 };
    AltToDMS(radalt_dms, radcenter.alt);

    //
    //寻找雷达作用范围,因为经纬高格子数应该一样，使用高程计算出格子数即可
    //
    float minalt = radcenter.alt - range;
    float minalt_dms[3] = { 0 };//最小高程
    AltToDMS(minalt_dms, minalt);
    float maxalt = radcenter.alt + range;
    float maxalt_dms[3] = { 0 };//最大高程
    AltToDMS(maxalt_dms, maxalt);

    float minlon_de = radlon_de - (range / (110000 * cos(radlat_de / 180 * PI)));
    float minlon_dms[3] = { 0 };//最小经度
    DecimalToDMS(minlon_dms, minlon_de);

    float minlat_de = radlat_de - (range / 110000);
    float minlat_dms[3] = { 0 };//最小纬度
    DecimalToDMS(minlat_dms, minlat_de);


    //确定移动步长
    auto dstep_dms = LevelStepMap[level];
    auto dstep_d = dstep_dms[0];//以度作为步长

    int edge_n = (maxalt_dms[0] - minalt_dms[0]) / dstep_d;

    int grid_n = edge_n * edge_n * edge_n;  //为什么三个维度上的网格数是一样的?
    Position* points = (Position*)malloc(sizeof(Position) * grid_n);
    if (points == NULL)
    {
        std::cout << "内存分配不成功!" << std::endl;
        exit(-1);
    }
    float alt[3] = { 0 };
    float lat[3] = { 0 };
    float lon[3] = { 0 };
    for (int i = 0; i < edge_n; i++) //高程
    {
        alt[0] = minalt_dms[0] + dstep_d * i;
        alt[1] = 0;       //其实不用每一步都赋值，只是为了后面更通用
        alt[2] = 0;
        for (int j = 0; j < edge_n; j++) //纬度
        {

            lat[0] = minlat_dms[0] + dstep_d * j;
            lat[1] = 0;
            lat[2] = 0;
            for (int z = 0; z < edge_n; z++) //经度
            {
                lon[0] = minlon_dms[0] + dstep_d * z;
                lon[1] = 0;
                lon[2] = 0;
                auto& pos = points[i * edge_n * edge_n + j * edge_n + z];
                pos.alt = DMSToAlt(alt);
                pos.lat[0] = lat[0];
                pos.lat[1] = lat[1];
                pos.lat[2] = lat[2];
                pos.lon[0] = lon[0];
                pos.lon[1] = lon[1];
                pos.lon[2] = lon[2];
            }
        }
    }
    return std::pair<Position*, int>{points, grid_n};
}
//10到15级，以分作步长
inline std::pair<Position*, int> CreateGridbyMi(Position radcenter, float range, unsigned short level)
{
    float radlat_de = DMSToDecimal(radcenter.lat);
    float radlon_de = DMSToDecimal(radcenter.lon);
    //先将radcenter的高程转化为度分秒的形式
    float radalt_dms[3] = { 0 };
    AltToDMS(radalt_dms, radcenter.alt);

    //
    //寻找雷达作用范围
    //
    float minalt = radcenter.alt - range;
    float minalt_dms[3] = { 0 };//最小高程
    AltToDMS(minalt_dms, minalt);
    float maxalt = radcenter.alt + range;
    float maxalt_dms[3] = { 0 };//最大高程
    AltToDMS(maxalt_dms, maxalt);

    float minlon_de = radlon_de - (range / (110000 * cos(radlat_de / 180 * PI)));
    float minlon_dms[3] = { 0 };//最小经度
    DecimalToDMS(minlon_dms, minlon_de);


    float minlat_de = radlat_de - (range / 110000);
    float minlat_dms[3] = { 0 };//最小纬度
    DecimalToDMS(minlat_dms, minlat_de);


    //确定移动步长
    auto dstep_dms = LevelStepMap[level];
    int dstep_m = dstep_dms[1];//以分作为步长（改成int，因为dms都应该是int型）

    //使用三个值：quot——商、remd——余数、flag——余数是否为0
    //计算出1度包含多少个步长quot以及是否有余数remd，flag用于标记能否整除。
    int quot = 60 / dstep_m;
    int remd = 60 % dstep_m;
    int flag = (remd == 0) ? 0 : 1;//余数为0，即能整除时标记为0，不能整除标记为1

    //先将min、max规范化到对应的网格
    minalt_dms[1] = (int(minalt_dms[1]) / dstep_m) * dstep_m;
    minlon_dms[1] = (int(minlon_dms[1]) / dstep_m) * dstep_m;
    minlat_dms[1] = (int(minlat_dms[1]) / dstep_m) * dstep_m;
    maxalt_dms[1] = (int(maxalt_dms[1]) / dstep_m) * dstep_m;


    // //计算网格的数量，以分为单位时，等于：度差×每度的格子数+（min对应的度内的起始格子数）+（max对应的度内的终止格子数）
    // int edge_n = (maxalt_dms[0] - minalt_dms[0]) * (quot + flag) + (quot + flag - int(minalt_dms[1]) / dstep_m) + (int(maxalt_dms[1]) / dstep_m);

    int realde2mi = (quot + flag) * dstep_m;//按照步长一度所对应的忽略边界问题的分
    //全部转化成分来计算(因为已经规范化了 所以直接相减就ok了)换算时一度对应的应该是按照步长一度所对应忽略边界问题的分
    int edge_n = int((maxalt_dms[0] - minalt_dms[0]) * realde2mi + maxalt_dms[1] - minalt_dms[1]) / dstep_m;

    int grid_n = edge_n * edge_n * edge_n;
    Position* points = (Position*)malloc(sizeof(Position) * grid_n);
    if (points == NULL)
    {
        std::cout << "内存分配不成功!" << std::endl;
        exit(-1);
    }

    //建立grids
    float alt[3] = { 0 };
    float lat[3] = { 0 };
    float lon[3] = { 0 };
    
    for (int i = 0; i < edge_n; i++) //高程
    {
        //解决边界问题
        int alt_plus = i * dstep_m;
        int alt_plus_d = (alt_plus + int(minalt_dms[1])) / realde2mi;
        int alt_plus_m = (alt_plus + int(minalt_dms[1])) % realde2mi;

        alt[0] = minalt_dms[0] + alt_plus_d;
        alt[1] = alt_plus_m;
        alt[2] = 0;
        for (int j = 0; j < edge_n; j++) //纬度
        {
            //解决边界问题
            int lat_plus = j * dstep_m;                      //对应示例：   8   16  24  32  40  48  56  64 
            int lat_plus_d = (lat_plus + int(minlat_dms[1])) / realde2mi;//度:0   0   0   0   0   0   0   1
            int lat_plus_m = (lat_plus + int(minlat_dms[1])) % realde2mi;//分:8   16  24  32  40  48  56  0

            lat[0] = minlat_dms[0] + lat_plus_d;
            lat[1] = lat_plus_m;
            lat[2] = 0;
            for (int z = 0; z < edge_n; z++) //经度
            {
                //解决边界问题
                int lon_plus = z * dstep_m;
                int lon_plus_d = (lon_plus + int(minlat_dms[1])) / realde2mi;
                int lon_plus_m = (lon_plus + int(minlat_dms[1])) % realde2mi;

                lon[0] = minlon_dms[0] + lon_plus_d;
                lon[1] = lon_plus_m;
                lon[2] = 0;
                auto& pos = points[i * edge_n * edge_n + j * edge_n + z];
                pos.alt = DMSToAlt(alt);
                pos.lat[0] = lat[0];
                pos.lat[1] = lat[1];
                pos.lat[2] = lat[2];
                pos.lon[0] = lon[0];
                pos.lon[1] = lon[1];
                pos.lon[2] = lon[2];
            }
        }
    }
    return std::pair<Position*, int>{points, grid_n};
}
inline std::pair<Position*, int> CreateGridbySe(Position radcenter, float range, unsigned short level)
{
     // --- 1. 计算探测范围的地理边界 (与原版类似) ---
     float radlat_de = DMSToDecimal(radcenter.lat);
     float radlon_de = DMSToDecimal(radcenter.lon);
 
     // 经纬度边界计算
     float minlon_de = radlon_de - (range / (110000 * cos(radlat_de / 180 * PI)));
     float maxlon_de = radlon_de + (range / (110000 * cos(radlat_de / 180 * PI)));
     float minlat_de = radlat_de - (range / 110000);
     float maxlat_de = radlat_de + (range / 110000);
 
     float minlon_dms[3] = {0}; DecimalToDMS(minlon_dms, minlon_de);
     float maxlon_dms[3] = {0}; DecimalToDMS(maxlon_dms, maxlon_de);
     float minlat_dms[3] = {0}; DecimalToDMS(minlat_dms, minlat_de);
     float maxlat_dms[3] = {0}; DecimalToDMS(maxlat_dms, maxlat_de);
 
     // --- ‼️ 核心修改 #1: 只计算雷达高度以上的区域 ---
     // 高程的最小值现在是雷达自身的高度，而不是 radcenter.alt - range
     float minalt = radcenter.alt;
     float maxalt = radcenter.alt + range;
     float minalt_dms[3] = {0}; AltToDMS(minalt_dms, minalt);
     float maxalt_dms[3] = {0}; AltToDMS(maxalt_dms, maxalt);
 
 
     // --- 2. 确定移动步长和修正参数 (与原版相同) ---
     auto dstep_dms = LevelStepMap[level];
     int dstep_s = dstep_dms[2]; // 以秒作为步长
 
     int quot = 60 / dstep_s;
     int remd = 60 % dstep_s;
     int flag = (remd == 0) ? 0 : 1;
     int realde2mi = (quot + flag) * dstep_s; // 一分对应的虚拟秒数
 
 
     // --- 3. ‼️ 核心修改 #2: 为每个轴向分别计算步数 (修复原版逻辑缺陷) ---
    // --- 3. 核心修正：将您的步数计算逻辑封装成一个可重用的函数 ---
    auto calculate_steps = [&](const float min_dms[3], const float max_dms[3]) -> int {
        // 将度、分、秒都转换为整数，避免浮点数精度问题
        int min_d = static_cast<int>(min_dms[0]);
        int min_m = static_cast<int>(min_dms[1]);
        int min_s = static_cast<int>(round(min_dms[2]));
        int max_d = static_cast<int>(max_dms[0]);
        int max_m = static_cast<int>(max_dms[1]);
        int max_s = static_cast<int>(round(max_dms[2]));

        // 将总范围转换为"虚拟秒"的总数
        int total_minutes = (max_d - min_d) * 60 + (max_m - min_m);
        long long total_virtual_seconds = total_minutes * realde2mi + (max_s - min_s);
        
        // 用总的虚拟秒数除以步长，向上取整
        return static_cast<int>(ceil(static_cast<double>(total_virtual_seconds) / dstep_s));
    };
 
     int steps_alt = calculate_steps(minalt_dms, maxalt_dms);
     int steps_lat = calculate_steps(minlat_dms, maxlat_dms);
     int steps_lon = calculate_steps(minlon_dms, maxlon_dms);
     
     // 如果任何一个方向步数为0或负数，则直接返回空
     if (steps_alt <= 0 || steps_lat <= 0 || steps_lon <= 0) {
         return {};
     }
         // 计算总点数，防止整数溢出
    size_t total_points = 0;
    if (steps_alt > 0 && steps_lat > 0 && steps_lon > 0) {
        // 先计算两维乘积，再乘第三维，避免溢出
        size_t points_per_alt = static_cast<size_t>(steps_lat) * steps_lon;
        total_points = static_cast<size_t>(steps_alt) * points_per_alt;
        
        // 检查溢出
        if (points_per_alt > 0 && total_points / points_per_alt != steps_alt) {
            std::cerr << "Warning: Integer overflow in grid points calculation" << std::endl;
            return {nullptr, 0}; // 溢出时返回空
        }
    }
    // 动态分配内存，而不是使用局部vector
    Position* result_points = nullptr;
    try {
        result_points = new Position[total_points];
    } catch (const std::bad_alloc& e) {
        std::cerr << "Error allocating memory for grid points: " << e.what() << std::endl;
        return {nullptr, 0};
    }
    size_t index = 0;
     float alt[3] = {0};
     float lat[3] = {0};
     float lon[3] = {0};
     for (int i = 0; i < steps_alt; i++) {
        int alt_plus = i * dstep_s;
        int alt_plus_m = (alt_plus + static_cast<int>(minalt_dms[2])) / realde2mi;
        int alt_plus_d = (alt_plus_m + static_cast<int>(minalt_dms[1])) / 60;
        alt[0] = minalt_dms[0] + alt_plus_d;
        alt[1] = (static_cast<int>(minalt_dms[1]) + alt_plus_m) % 60;
        alt[2] = (alt_plus + static_cast<int>(minalt_dms[2])) % realde2mi;

        for (int j = 0; j < steps_lat; j++) {
            // ... (与 alt 类似的逻辑计算 lat[3]) ...
            int lat_plus = j * dstep_s;
            int lat_plus_m = (lat_plus + static_cast<int>(minlat_dms[2])) / realde2mi;
            int lat_plus_d = (lat_plus_m + static_cast<int>(minlat_dms[1])) / 60;
            lat[0] = minlat_dms[0] + lat_plus_d;
            lat[1] = (static_cast<int>(minlat_dms[1]) + lat_plus_m) % 60;
            lat[2] = (lat_plus + static_cast<int>(minlat_dms[2])) % realde2mi;

            for (int z = 0; z < steps_lon; z++) {
                // ... (与 alt 类似的逻辑计算 lon[3]) ...
                int lon_plus = z * dstep_s;
                int lon_plus_m = (lon_plus + static_cast<int>(minlon_dms[2])) / realde2mi;
                int lon_plus_d = (lon_plus_m + static_cast<int>(minlon_dms[1])) / 60;
                lon[0] = minlon_dms[0] + lon_plus_d;
                lon[1] = (static_cast<int>(minlon_dms[1]) + lon_plus_m) % 60;
                lon[2] = (lon_plus + static_cast<int>(minlon_dms[2])) % realde2mi;

                Position pos;
                pos.alt = DMSToAlt(alt);
                // ... 赋值 pos.lat 和 pos.lon ...
                result_points[index].alt = DMSToAlt(alt);
                memcpy(result_points[index].lat, lat, sizeof(float) * 3);
                memcpy(result_points[index].lon, lon, sizeof(float) * 3);
                index++;
            }
        }
    }
    return {result_points, static_cast<int>(total_points)};
}

// 在 grid.hpp 中，添加到 CreateGridbySe 函数之后// 1. 先创建一个可重用的辅助函数，用于将度分秒转换为总秒数
inline double DMSToTotalSeconds(const float dms[3]) {
    double sign = (dms[0] < 0 || dms[1] < 0 || dms[2] < 0) ? -1.0 : 1.0;
    return sign * (std::abs(dms[0]) * 3600.0 + std::abs(dms[1]) * 60.0 + std::abs(dms[2]));
}
// 2. 再创建一个辅助函数，用于将总秒数转换回度分秒
inline void TotalSecondsToDMS(float dms[3], double total_seconds) {
    double sign = (total_seconds < 0) ? -1.0 : 1.0;
    total_seconds = std::abs(total_seconds);

    dms[0] = static_cast<int>(total_seconds / 3600.0);
    total_seconds -= dms[0] * 3600.0;
    dms[1] = static_cast<int>(total_seconds / 60.0);
    dms[2] = total_seconds - dms[1] * 60.0;
    
    dms[0] *= sign;
    // 分和秒通常保持正值，符号由“度”体现，但为保持一致性，此处也应用符号// 在实际使用中可能需要根据具体需求调整符号处理逻辑
    dms[1] *= sign;
    dms[2] *= sign;
}
// 新增的函数，用于处理22级及以上的高精度剖分
inline std::pair<Position*, int> CreateGridByFractionalSec(Position radcenter, float range, unsigned short level){
    // --- 1. 动态计算步长 (单位：秒) ---// level 22 -> 1/2s, level 23 -> 1/4s, etc.
    double step_s = 1.0 / (1 << (level - 21));

    // --- 2. 计算探测范围的地理边界 (与CreateGridbySe逻辑相同) ---
    float radlat_de = DMSToDecimal(radcenter.lat);
    float radlon_de = DMSToDecimal(radcenter.lon);

    float minlon_de = radlon_de - (range / (110000.0f * cos(radlat_de / 180.0f * PI)));
    float maxlon_de = radlon_de + (range / (110000.0f * cos(radlat_de / 180.0f * PI)));
    float minlat_de = radlat_de - (range / 110000.0f);
    float maxlat_de = radlat_de + (range / 110000.0f);
    float minalt = radcenter.alt - 50; //TODO 之后可能要改正，因为只能探测雷达高度以上的区域
    float maxalt = radcenter.alt + range;

    float minlon_dms[3] = {0}; DecimalToDMS(minlon_dms, minlon_de);
    float maxlon_dms[3] = {0}; DecimalToDMS(maxlon_dms, maxlon_de);
    float minlat_dms[3] = {0}; DecimalToDMS(minlat_dms, minlat_de);
    float maxlat_dms[3] = {0}; DecimalToDMS(maxlat_dms, maxlat_de);
    float minalt_dms[3] = {0}; AltToDMS(minalt_dms, minalt);
    float maxalt_dms[3] = {0}; AltToDMS(maxalt_dms, maxalt);

    // --- 3. 计算每个轴向的总步数 ---
    auto calculate_steps = [&](const float min_dms[3], const float max_dms[3]) -> int {
        double min_total_sec = DMSToTotalSeconds(min_dms);
        double max_total_sec = DMSToTotalSeconds(max_dms);
        double total_range_sec = max_total_sec - min_total_sec;
        if (total_range_sec <= 0) return 0;
        return static_cast<int>(ceil(total_range_sec / step_s)) + 1;
    };

    int steps_alt = calculate_steps(minalt_dms, maxalt_dms);
    int steps_lat = calculate_steps(minlat_dms, maxlat_dms);
    int steps_lon = calculate_steps(minlon_dms, maxlon_dms);

    if (steps_alt <= 0 || steps_lat <= 0 || steps_lon <= 0) {
      return {nullptr, 0};
    }

    // --- 4. 分配内存并进行溢出检查 ---
    size_t total_points = static_cast<size_t>(steps_alt) * static_cast<size_t>(steps_lat) * static_cast<size_t>(steps_lon);
    std::cout << "Debug: total_points to be allocated: " << total_points << std::endl;
    if (total_points > 200000000) { // 设置一个上限，防止内存爆炸std::cerr << "Error: Grid size (" << total_points << ") is too large. Aborting." << std::endl;
        return {nullptr, 0};
    }
    
    Position* points = new Position[total_points];

    // --- 5. 循环生成网格点 ---
    double min_alt_total_sec = DMSToTotalSeconds(minalt_dms);
    double min_lat_total_sec = DMSToTotalSeconds(minlat_dms);
    double min_lon_total_sec = DMSToTotalSeconds(minlon_dms);
    
    size_t index = 0;
    for (int i = 0; i < steps_alt; ++i) {
        float current_alt_dms[3];
        TotalSecondsToDMS(current_alt_dms, min_alt_total_sec + i * step_s);

        for (int j = 0; j < steps_lat; ++j) {
            float current_lat_dms[3];
            TotalSecondsToDMS(current_lat_dms, min_lat_total_sec + j * step_s);

            for (int k = 0; k < steps_lon; ++k) {
                float current_lon_dms[3];
                TotalSecondsToDMS(current_lon_dms, min_lon_total_sec + k * step_s);

                auto& pos = points[index++];
                pos.alt = DMSToAlt(current_alt_dms);
                memcpy(pos.lat, current_lat_dms, sizeof(float) * 3);
                memcpy(pos.lon, current_lon_dms, sizeof(float) * 3);
            }
        }
    }

    return {points, static_cast<int>(total_points)};
}




//16到21级，以秒作步长
// std::pair<Position*, int> CreateGridbySe(Position radcenter, float range, unsigned short level)
// {
//     float radlat_de = DMSToDecimal(radcenter.lat);
//     float radlon_de = DMSToDecimal(radcenter.lon);
//     //先将radcenter的高程转化为度分秒的形式
//     float radalt_dms[3] = { 0 };
//     AltToDMS(radalt_dms, radcenter.alt);

//     //
//     //寻找雷达作用范围
//     //
//     float minalt = radcenter.alt - range;
//     float minalt_dms[3] = { 0 };//最小高程
//     AltToDMS(minalt_dms, minalt);
//     float maxalt = radcenter.alt + range;
//     float maxalt_dms[3] = { 0 };//最大高程
//     AltToDMS(maxalt_dms, maxalt);

//     float minlon_de = radlon_de - (range / (110000 * cos(radlat_de / 180 * PI)));
//     float minlon_dms[3] = { 0 };//最小经度
//     DecimalToDMS(minlon_dms, minlon_de);


//     float minlat_de = radlat_de - (range / 110000);
//     float minlat_dms[3] = { 0 };//最小纬度
//     DecimalToDMS(minlat_dms, minlat_de);

//     //确定移动步长
//     auto dstep_dms = LevelStepMap[level];
//     int dstep_s = dstep_dms[2];//以秒作为步长（改成int，因为dms都应该是int型）

//     //使用三个值：quot——商、remd——余数、flag——余数是否为0
//     int quot = 60 / dstep_s;
//     int remd = 60 % dstep_s;
//     int flag = (remd == 0) ? 0 : 1;//余数为0，即能整除时标记为0，不能整除标记为1

//     //先将min、max规范化到对应的网格
//     minalt_dms[2] = (int(minalt_dms[2]) / dstep_s) * dstep_s;
//     minlon_dms[2] = (int(minlon_dms[2]) / dstep_s) * dstep_s;
//     minlat_dms[2] = (int(minlat_dms[2]) / dstep_s) * dstep_s;
//     maxalt_dms[2] = (int(maxalt_dms[2]) / dstep_s) * dstep_s;

//     // 计算网格的数量，以秒为单位时，等于：分差×每分的格子数+（min对应的分内的起始格子数）+（max对应的分内的终止格子数）
//     // int edge_n = ((maxalt_dms[0] * 60 + maxalt_dms[1]) - (minalt_dms[0] * 60 + minalt_dms[1])) * (quot + flag) + (quot + flag - int(minalt_dms[2]) / dstep_s) + (int(maxalt_dms[2]) / dstep_s);
    
    
//     int realde2mi = (quot + flag) * dstep_s;//按照步长一分所对应的忽略边界问题的秒

//     //全部转化成分来计算(因为已经规范化了 所以直接相减就ok了
//     int edge_n = int((maxalt_dms[0] - minalt_dms[0]) *60*realde2mi  + (maxalt_dms[1] - minalt_dms[1])*realde2mi+maxalt_dms[0]-minalt_dms[0]) / dstep_s;

//     int grid_n = edge_n * edge_n * edge_n;
//     Position* points = (Position*)malloc(sizeof(Position) * grid_n);
//     if (points == NULL)
//     {
//         std::cout << "内存分配不成功!" << std::endl;
//         exit(-1);
//     }

//     //建立grids
//     float alt[3] = { 0 };
//     float lat[3] = { 0 };
//     float lon[3] = { 0 };

//     for (int i = 0; i < edge_n; i++) //高程
//     {   //解决边界问题 
//         int alt_plus = i * dstep_s;
//         int alt_plus_m = (alt_plus + int(minalt_dms[2])) / realde2mi;
//         int alt_plus_d = (alt_plus_m + int(minalt_dms[1])) / 60;
//         int alt_plus_s = (alt_plus + int(minalt_dms[2])) % realde2mi;

//         alt[0] = minalt_dms[0] + alt_plus_d;
//         alt[1] = (int(minalt_dms[1]) + alt_plus_m) % 60;
//         alt[2] = alt_plus_s;
//         for (int j = 0; j < edge_n; j++) //纬度
//         {
//             //解决边界问题 
//             int lat_plus = j * dstep_s;
//             int lat_plus_m = (lat_plus + int(minlat_dms[2])) / realde2mi;
//             int lat_plus_d = (lat_plus_m + int(minlat_dms[1])) / 60;
//             int lat_plus_s = (lat_plus + int(minlat_dms[2])) % realde2mi;

//             lat[0] = minlat_dms[0] + lat_plus_d;
//             lat[1] = (int(minlat_dms[1]) + lat_plus_m) % 60;
//             lat[2] = lat_plus_s;
//             for (int z = 0; z < edge_n; z++) //经度
//             {
//                 //解决边界问题 
//                 int lon_plus = z * dstep_s;
//                 int lon_plus_m = (lon_plus + int(minlon_dms[2])) / realde2mi;
//                 int lon_plus_d = (lon_plus_m + int(minlon_dms[1])) / 60;
//                 int lon_plus_s = (lon_plus + int(minlon_dms[2])) % realde2mi;

//                 lon[0] = minlon_dms[0] + lon_plus_d;
//                 lon[1] = (int(minlon_dms[1]) + lon_plus_m) % 60;
//                 lon[2] = lon_plus_s;
//                 auto& pos = points[i * edge_n * edge_n + j * edge_n + z];
//                 pos.alt = DMSToAlt(alt);
//                 pos.lat[0] = lat[0];
//                 pos.lat[1] = lat[1];
//                 pos.lat[2] = lat[2];
//                 pos.lon[0] = lon[0];
//                 pos.lon[1] = lon[1];
//                 pos.lon[2] = lon[2];
//             }
//         }
//     }
//     return std::pair<Position*, int>{points, grid_n};
// }
//计算网格
//input:
//	radcenter:雷达中心点
//	range:雷达作用距离,单位为米
// .level:网格层级
//output:
//	Position[1000000]：遍历得到的100万个点的经纬高

//不同单位会增加理解的难度，因此将高程转化到度分秒，高度上，地心的高程为0度，512度对应了512×110km即56302米
inline std::pair<Position *, int> CreateGrid(Position radcenter, float range, unsigned short level)
{
    std::pair<Position *,int> point_list{nullptr,0};
    if(level<10)//level<10时，步长单位为度
    {
        point_list = CreateGridbyDe(radcenter, range, level);
    }
    else if(level<16)
    {
        point_list = CreateGridbyMi(radcenter, range, level);
    }
    else if(level<22)
    {
        point_list = CreateGridbySe(radcenter, range, level);
    }else {
        point_list = CreateGridByFractionalSec(radcenter, range, level);//注意最高到26级，不能超过26级，26级的时候网格的边长已经是1m了已经足够的精细
    }
    return point_list;
}
