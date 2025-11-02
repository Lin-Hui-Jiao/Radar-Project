#include "hiradar/radar.hpp"
#include <cmath>
#include <proj.h>
#include "hiradar/grid.hpp"
#include <png.h> // 用于PNG操作
#include <sstream> // 用于std::ostringstream
// 一个简单的结构体，用于表示RGB颜色
struct Color { unsigned char r, g, b; };
// 首先定义一个结构体来表示色带上的一个颜色点
struct ColorPoint {
    float value; // 归一化值 [0.0, 1.0]
    int r, g, b; // 对应的RGB颜色
};

// 【核心】定义雷达能量强度颜色查找表
// 设计理念：低能量=浅黄色，高能量=纯红色，突出能量差异对比
// 符合用户需求：最强能量=鲜红色，最弱能量=浅黄色，无信号=灰色
static const std::vector<ColorPoint> colormap = {
    {0.0f, 255, 255, 180}, // 0.0: 浅黄色 - 代表最低有效信号
    {0.1f, 255, 240, 160}, // 0.1: 淡黄色
    {0.2f, 255, 220, 140}, // 0.2: 黄色
    {0.3f, 255, 200, 100}, // 0.3: 橙黄色
    {0.4f, 255, 170,  60}, // 0.4: 橙色
    {0.5f, 255, 140,  40}, // 0.5: 深橙色
    {0.6f, 255, 110,  20}, // 0.6: 橙红色
    {0.7f, 255,  80,  10}, // 0.7: 红橙色
    {0.8f, 255,  50,   5}, // 0.8: 亮红色
    {0.9f, 255,  25,   2}, // 0.9: 深红色
    {1.0f, 255,   0,   0}  // 1.0: 纯红色 - 代表最高信号强度
};

Color lerp(const ColorPoint& c1, const ColorPoint& c2, float value) {
    float t = (value - c1.value) / (c2.value - c1.value);
    unsigned char r = static_cast<unsigned char>(c1.r + t * (c2.r - c1.r));
    unsigned char g = static_cast<unsigned char>(c1.g + t * (c2.g - c1.g));
    unsigned char b = static_cast<unsigned char>(c1.b + t * (c2.b - c1.b));
    return {r, g, b};
}


Color mapDensityToColorLog_v2(float density, float min_density, float max_density) {
    // --- 1. 特殊情况：能量为0或负值时显示灰色 ---
    if (density <= 0) {
        return {128, 128, 128}; // 灰色表示无信号区域
    }

    // --- 2. 参数有效性检查 ---
    if (max_density <= min_density || max_density <= 0 || min_density <= 0) {
        return {128, 128, 128}; // 参数无效时返回灰色
    }

    // --- 3. 对数映射：将线性功率密度转换为对数尺度 ---
    // 使用自然对数而不是分贝，获得更平滑的过渡
    float log_density = std::log(density);
    float log_min = std::log(min_density);
    float log_max = std::log(max_density);

    // 检查对数范围有效性
    if (std::abs(log_max - log_min) < 1e-6) {
        // 如果最大值和最小值相等，返回中等强度的颜色
        return {static_cast<unsigned char>(colormap[colormap.size()/2].r),
                static_cast<unsigned char>(colormap[colormap.size()/2].g),
                static_cast<unsigned char>(colormap[colormap.size()/2].b)};
    }

    // --- 4. 归一化到 [0, 1] 区间 ---
    float normalized_v = (log_density - log_min) / (log_max - log_min);

    // 限制范围，确保在 [0, 1] 内
    normalized_v = std::max(0.0f, std::min(1.0f, normalized_v));

    // --- 5. 应用幂函数调整，强化高能量区域的对比度 ---
    // 使用平方根函数，使低能量区域颜色变化更细腻
    normalized_v = std::sqrt(normalized_v);

    // --- 6. 在颜色表中查找并插值 ---
    // 处理边界情况
    if (normalized_v <= colormap.front().value) {
        return {static_cast<unsigned char>(colormap.front().r),
                static_cast<unsigned char>(colormap.front().g),
                static_cast<unsigned char>(colormap.front().b)};
    }
    if (normalized_v >= colormap.back().value) {
        return {static_cast<unsigned char>(colormap.back().r),
                static_cast<unsigned char>(colormap.back().g),
                static_cast<unsigned char>(colormap.back().b)};
    }

    // 找到合适的颜色区间并进行插值
    for (size_t i = 0; i < colormap.size() - 1; ++i) {
        const auto& c1 = colormap[i];
        const auto& c2 = colormap[i + 1];
        if (normalized_v >= c1.value && normalized_v <= c2.value) {
            return lerp(c1, c2, normalized_v);
        }
    }

    // 备用返回值（理论上不会到达）
    return {static_cast<unsigned char>(colormap.back().r),
            static_cast<unsigned char>(colormap.back().g),
            static_cast<unsigned char>(colormap.back().b)};
}
/**
 * 这是一个底层的单点计算函数。它的唯一职责是，
 * 根据一个点相对于雷达天线朝向的相对球坐标 {theta, phi, R}，计算出这个点的电磁功率密度。
 */
float Radar::PointPowerDensity(float tpr[3]) // tpr[3]这是球坐标系
{
    caltools ct;
    if (tpr[2] < 0) //R为负就不进行计算
    {
        return 0;
    }
    // a. 根据雷达频率(_params.f)和信号俯仰角(tpr[1])，查表并插值计算出大气衰减的分贝值(dB/km)
    float reduce_DB = ct.choose_rDB(_params.f, tpr[1]); 
    // b. 将dB/km单位的衰减值，转换为一个线性的衰减系数 alpha
    float alpha = ct.cal_alpha(reduce_DB);
    // c. 根据衰减系数alpha和距离(tpr[2])，计算出总的大气损耗倍数 La
    float La = ct.cal_La(alpha, tpr[2]);
    // 查表得到方位角的方向系数
    float etheta = _etheta_table[int(tpr[0] * PRECISION)];
    // 查表得到俯仰角的方向系数
    float ephi = _ephi_table[int(tpr[1] * PRECISION)];
    float Gthetaphi = ct.cal_Gthetaphi(_params.G, etheta, ephi);
    float qt = ct.cal_Qt(_params.Pt, tpr[2], La, Gthetaphi);
    return qt;
}
/**
 * 计算在雷达当前固定的姿态和朝向下，一系列网格点（pos_list）的功率密度分布。
 * 可以理解为对雷达某一瞬间的“快照”进行计算。
 */
//这个函数之后还是需要进行改正
void Radar::PowerDensity(const Position *pos_list, size_t n)
{
//     int p_n, p_id;
//     MPI_Comm_size(MPI_COMM_WORLD, &p_n);
//     MPI_Comm_rank(MPI_COMM_WORLD, &p_id);
//     auto chunk = n / p_n;
//     auto start = p_id * chunk;
//     auto end = p_id < p_n - 1 ? (p_id + 1) * chunk : n;
//     auto local_n = end - start;
//     float *values = (float *)malloc(sizeof(float) * n);
//     caltools ct;
//     float(*tpr)[3] = (float(*)[3])malloc(sizeof(float) * n * 3);
//     time_t start_time = clock();
//     //坐标变换流水线
//     for (int i = start; i < end; i++)
//     {
//         float xyz[3]{};
//         float thephiR[3]{};
//         //全局坐标系转为雷达所在位置为原点，正北方向为x轴
//         ct.lon2xyz(_pos.lon, _pos.lat, _pos.alt, pos_list[i].lon, pos_list[i].lat, pos_list[i].alt, xyz);
//         //变换坐标系之后，点在新坐标系中的位置(x,y,z)
//         ct.trans_axes(xyz, _base_theta + _theta, _base_phi + _phi);
//         // 笛卡尔坐标系转球坐标系
//         ct.xyz2tpr(xyz, thephiR);
//         //现在thephiR是p点与载具加上雷达方向的夹角
//         float radxR = thephiR[2] / 1000; //米转换为千米

//         tpr[i][0] = thephiR[0];
//         tpr[i][1] = thephiR[1];
//         tpr[i][2] = radxR;
//     }

//     //这里是错误
//     for (int i = 0; i < n; i++)
//     {
//         values[i] = PointPowerDensity(tpr[i]);
//     }

//     time_t end_time = clock();
// #ifdef DEBUG
//     if (p_id == 0)
//     {
//         cout << "网格数：" << n << "，计算用时：";
//         cout << (end_time - start_time) << endl;
//     }
// #endif

//     free(tpr);
//     _writer->Write(pos_list + start, values + start, local_n, _timestamp);//并发的写入到数据库中
//     free(values);
}
/**
 * 它计算的不是雷达某一瞬间的功率密度，而是当雷达天线在其允许的扫描范围内扫视时，
 * 每个网格点可能接收到的最大功率密度。它生成的是一张“最优探测能力图”。
 */
void Radar::CapablePowerDensity(float out[], const Position *pos_list, size_t n, float radius[])
{
   
    #pragma omp parallel
    {
        PROJ_Manager proj_manager;
        caltools ct;
        #pragma omp for schedule(dynamic)
        for(size_t i = 0; i < n; i++){
            const auto& target_pos = pos_list[i];
            // a. 遮挡判断
            if (isOccluded(target_pos,&proj_manager))
            {
                // 如果被遮挡，直接将该点的输出能量置零
                out[i] = 0.0f;
                if (radius != NULL) {
                    radius[i] = -1.0f; // 标记距离为无效
                }
                continue; // 使用 continue 跳过当前循环，处理下一个点
            }
            // b. 坐标转换和能量计算 (这部分逻辑完全不变)
            float xyz[3]{};
            float thephiR[3]{};
            ct.lon2xyz(_pos.lon, _pos.lat, _pos.alt, target_pos.lon, target_pos.lat, target_pos.alt, xyz);
            ct.trans_axes(xyz, _base_theta, _base_phi);
            ct.xyz2tpr(xyz, thephiR);
            float minradxtheta = 0;
            float minradxphi = 0;
            float radxR = thephiR[2] / 1000.0f; 

            if (radius != NULL) {
                radius[i] = radxR;
            }
            //先计算每点与载具角度的差值，再看其是否在旋转范围内
            float radxtheta = thephiR[0];
            float radxphi = thephiR[1];
            radxtheta = ct.normalize_angle(radxtheta);
            //找到最小的theta
            if (radxtheta < _params.mintheta)
            {
                minradxtheta = _params.mintheta - radxtheta;
            }
            else if ((radxtheta >= _params.mintheta) && (radxtheta <= _params.maxtheta))
            {
                minradxtheta = 0;
            }
            else if (radxtheta > _params.maxtheta)
            {
                minradxtheta = radxtheta - _params.maxtheta;
            }
            //找到最小的phi
            if (radxphi < _params.minphi)
            {
                minradxphi = _params.minphi - radxphi;
            }
            else if ((radxphi >= _params.minphi) && (radxphi <= _params.maxphi))
            {
                minradxphi = 0;
            }
            else if (radxphi > _params.maxphi)
            {
                minradxphi = radxphi - _params.maxphi;

            }
            
            float tpr[3];
            if (minradxtheta >= 90 || minradxphi >= 90)
            {
                tpr[0] = 0;
                tpr[1] = 0;
                tpr[2] = -1;
            }
            else
            {
                tpr[0] = minradxtheta;
                tpr[1] = minradxphi;
                tpr[2] = radxR;
            }   
            out[i] = PointPowerDensity(tpr);
        }
    }
}
void Radar::CapableGroundPowerDensity(int range,double step){
    std::cout << "--> Planning 2D ground grid..." << std::endl;
    // a. 确定计算区域 (Area of Interest) 的局部坐标范围
     // a. 确定计算区域 (Area of Interest) 的局部坐标范围
     double min_x = _cached_local_x - range; //_cached_local_x和_cached_local_y是转换后的在香港坐标系下的雷达坐标
     double max_x = _cached_local_x + range;
     double min_y = _cached_local_y - range;
     double max_y = _cached_local_y + range;
     
     // b. 根据层级确定二维网格的步长 (分辨率)
     // 这是一个简化的方式，您可以根据需要调整来控制网格密度
    
     
     int steps_x = static_cast<int>(ceil((max_x - min_x) / step));
     int steps_y = static_cast<int>(ceil((max_y - min_y) / step));
     size_t n = static_cast<size_t>(steps_x) * steps_y;
 
     if (n == 0) {
         std::cout << "Error: Grid size is zero." << std::endl;
         return;
     }
     std::cout << "--> Grid planned with " << n << " points." << std::endl;
 
     // --- 2. 分配内存 ---
     Position* pos_list = new Position[n];      // 存储最终包含海拔的3D坐标
     float* density_values = new float[n]; // 存储计算结果
 
     // --- 3. OpenMP 并行计算 ---
     std::cout << "--> Starting parallel computation (including ground height query)..." << std::endl;
     time_t start_time = clock();
 
     #pragma omp parallel
     {
         // 每个线程拥有自己独立的PROJ和caltools实例
         PROJ_Manager proj_manager;
         caltools ct;
         PJ* proj_to_wgs84 = proj_create_crs_to_crs(proj_manager.ctx, "EPSG:2326", "EPSG:4326", NULL);
         #pragma omp for schedule(dynamic)
         for (int i = 0; i < steps_y; ++i) {
             for (int j = 0; j < steps_x; ++j) {
                 size_t index = static_cast<size_t>(i) * steps_x + j;
 
                 // a. 计算当前点的局部投影坐标 (x, y)
                 double target_local_x = min_x + j * step;
                 double target_local_y = min_y + i * step;
 
                 // b. 【关键】查询该点的真实地表高度 z
                 double target_terrain_alt_m = _rtree->Getheight3d(target_local_x - 1975, 43 + target_local_y, HeightCallback) + 0.1;
                 // c. 构造完整的 Position 对象，用于写入数据库
                 pos_list[index].alt = target_terrain_alt_m;
                 // 将局部坐标反算为经纬度
                 PJ_COORD local_coord = proj_coord(target_local_y + 815500, target_local_x + 835000, 0, 0);
                 PJ_COORD wgs84_coord = proj_trans(proj_to_wgs84, PJ_FWD, local_coord);
                 DecimalToDMS(pos_list[index].lon, wgs84_coord.xy.y);
                 DecimalToDMS(pos_list[index].lat, wgs84_coord.xy.x);
                 
                 // 计算起始和结束点
                 double indexRange_x = -1975;
                 double indexRange_y = 43;
                 double target_abs_alt_m = target_terrain_alt_m;
                 double start_point[3] = {indexRange_x + _cached_local_x, indexRange_y + _cached_local_y, _cached_abs_alt_m};
                 double end_point[3]   = {indexRange_x + target_local_x,  indexRange_y + target_local_y,  target_abs_alt_m};
                 Rect3d search_rect(target_local_x + indexRange_x,indexRange_y + target_local_y, target_terrain_alt_m, indexRange_x + _cached_local_x, indexRange_y + _cached_local_y, _cached_abs_alt_m); // 250*250 500*500

                //  if(_rtree->Intersect3d(search_rect.min,search_rect.max, IntersectCallback)){
                if(_rtree->Intersect3d(search_rect.min,search_rect.max, IntersectCallback)){
                    density_values[index] = 0.0f; // 0 代表被遮挡
                    continue;
                }
                // b. 坐标转换和能量计算 (这部分逻辑完全不变)
                float xyz[3]{};
                float thephiR[3]{};
                ct.lon2xyz(_pos.lon, _pos.lat, _pos.alt, pos_list[index].lon, pos_list[index].lat, pos_list[index].alt, xyz);
                ct.trans_axes(xyz, _base_theta, _base_phi);
                ct.xyz2tpr(xyz, thephiR);
                float minradxtheta = 0;
                float minradxphi = 0;
                float radxR = thephiR[2] / 1000.0f; 
                //先计算每点与载具角度的差值，再看其是否在旋转范围内
                float radxtheta = thephiR[0];
                float radxphi = thephiR[1];
                radxtheta = ct.normalize_angle(radxtheta);
                //找到最小的theta
                if (radxtheta < _params.mintheta)
                {
                    minradxtheta = _params.mintheta - radxtheta;
                }
                else if ((radxtheta >= _params.mintheta) && (radxtheta <= _params.maxtheta))
                {
                    minradxtheta = 0;
                }
                else if (radxtheta > _params.maxtheta)
                {
                    minradxtheta = radxtheta - _params.maxtheta;
                }
                //找到最小的phi
                if (radxphi < _params.minphi)
                {
                    minradxphi = _params.minphi - radxphi;
                }
                else if ((radxphi >= _params.minphi) && (radxphi <= _params.maxphi))
                {
                    minradxphi = 0;
                }
                else if (radxphi > _params.maxphi)
                {
                    minradxphi = radxphi - _params.maxphi;

                }
                
                float tpr[3];
                if (minradxtheta >= 90 || minradxphi >= 90)
                {
                    tpr[0] = 0;
                    tpr[1] = 0;
                    tpr[2] = -1;
                }
                else
                {
                    tpr[0] = minradxtheta;
                    tpr[1] = minradxphi;
                    tpr[2] = radxR;
                }   
                density_values[index] = PointPowerDensity(tpr);
                //density_values[index] = 1;
             }
         }
         proj_destroy(proj_to_wgs84);
     } // 并行区域结束
      // --- 4. 颜色映射 ---
    std::cout << "--> Mapping densities to colors..." << std::endl;
    float min_visible_density = 1e12, max_visible_density = 0.0f;
    for (size_t i = 0; i < n; ++i) {
        if (density_values[i] > 0) { // 只在可见点中寻找最大最小值
            if (density_values[i] < min_visible_density) min_visible_density = density_values[i];
            if (density_values[i] > max_visible_density) max_visible_density = density_values[i];
        }
    }
    if (min_visible_density < 1e-6f) { // 假设1e-6是一个合理的最小能量阈值
        min_visible_density = 1e-6f;
   }
    png_bytep* rowPointers = (png_bytep*)malloc(steps_y * sizeof(png_bytep));
    for (int i = 0; i < steps_y; i++)
        rowPointers[i] = (png_bytep)malloc(steps_x * 4 * sizeof(png_byte));
    for (int i = 0; i < steps_y; i++) {
        for (int j = 0; j < steps_x; j++) {
            size_t index = static_cast<size_t>(i) * steps_x + j;
            float density = density_values[index];
            Color color;
            color = mapDensityToColorLog_v2(density, min_visible_density, max_visible_density);
            int ii = steps_y - 1 - i; // 垂直翻转以匹配图像坐标系
            rowPointers[ii][j * 4 + 0] = color.r;
            rowPointers[ii][j * 4 + 1] = color.g;
            rowPointers[ii][j * 4 + 2] = color.b;
            rowPointers[ii][j * 4 + 3] = 255; // Alpha通道 (不透明)
        }
    }
    // --- 5. 保存PNG图像 ---
    std::string output_filename = "capable_ground_power_density.png";
    FILE *fp = fopen(output_filename.c_str(), "wb");
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    png_infop info_ptr = png_create_info_struct(png_ptr);
    png_init_io(png_ptr, fp);
    png_set_IHDR(png_ptr, info_ptr, steps_x, steps_y, 8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    png_write_info(png_ptr, info_ptr);
    png_write_image(png_ptr, rowPointers);
    png_write_end(png_ptr, NULL);
    fclose(fp);
    png_destroy_write_struct(&png_ptr, &info_ptr);

 
    // // 5. 将计算结果批量写入数据库
    // std::cout << "--> Writing " << n << " results to database..." << std::endl;
    // _writer->Write(pos_list, density_values, n, _timestamp);
    std::cout << "--> Write operation completed." << std::endl;
    for (int i = 0; i < steps_y; i++) {
        free(rowPointers[i]); // 释放每一行的内存
    }
    free(rowPointers); // 释放行指针数组本身的内存
    delete[] density_values;
    delete[] pos_list;
}

void Radar::CapableFixedHeightPowerDensity(int range, double step, double fixed_height) {
    std::cout << "--> Planning 2D fixed-height grid at " << fixed_height << "m altitude..." << std::endl;

    // 1. 高度合理性检查
    if (fixed_height < 1) {
        std::cerr << "Error: Fixed height cannot be negative (" << fixed_height << "m)" << std::endl;
        return;
    }

    // 确保指定高度不会太接近雷达高度，避免计算异常
    double radar_terrain_alt = _cached_abs_alt_m - _pos.alt; // 雷达下方地表高度
    if (fixed_height < radar_terrain_alt - 50) {
        std::cerr << "Warning: Fixed height (" << fixed_height << "m) is significantly below radar terrain ("
                  << radar_terrain_alt << "m). Results may be unreliable." << std::endl;
    }

    // 2. 确定计算区域的局部坐标范围（与CapableGroundPowerDensity相同）
    double min_x = _cached_local_x - range;
    double max_x = _cached_local_x + range;
    double min_y = _cached_local_y - range;
    double max_y = _cached_local_y + range;

    // 3. 根据步长计算网格步数
    int steps_x = static_cast<int>(ceil((max_x - min_x) / step));
    int steps_y = static_cast<int>(ceil((max_y - min_y) / step));
    size_t n = static_cast<size_t>(steps_x) * steps_y;

    if (n == 0) {
        std::cout << "Error: Grid size is zero." << std::endl;
        return;
    }
    std::cout << "--> Grid planned with " << n << " points at fixed height " << fixed_height << "m." << std::endl;

    // 4. 分配内存
    Position* pos_list = new Position[n];
    float* density_values = new float[n];

    // 5. OpenMP 并行计算
    std::cout << "--> Starting parallel computation for fixed-height grid..." << std::endl;
    time_t start_time = clock();

    #pragma omp parallel
    {
        // 每个线程拥有自己独立的PROJ和caltools实例
        PROJ_Manager proj_manager;
        caltools ct;
        PJ* proj_to_wgs84 = proj_create_crs_to_crs(proj_manager.ctx, "EPSG:2326", "EPSG:4326", NULL);

        #pragma omp for schedule(dynamic)
        for (int i = 0; i < steps_y; ++i) {
            for (int j = 0; j < steps_x; ++j) {
                size_t index = static_cast<size_t>(i) * steps_x + j;

                // a. 计算当前点的局部投影坐标
                double target_local_x = min_x + j * step;
                double target_local_y = min_y + i * step;

                // b. 【关键修改】使用固定高度而不是查询地表高度
                double target_fixed_alt_m = fixed_height;

                // c. 构造完整的 Position 对象
                pos_list[index].alt = target_fixed_alt_m;

                // 将局部坐标反算为经纬度
                PJ_COORD local_coord = proj_coord(target_local_y + 815500, target_local_x + 835000, 0, 0);
                PJ_COORD wgs84_coord = proj_trans(proj_to_wgs84, PJ_FWD, local_coord);
                DecimalToDMS(pos_list[index].lon, wgs84_coord.xy.y);
                DecimalToDMS(pos_list[index].lat, wgs84_coord.xy.x);

                // d. 遮挡检测（保持与原函数相同的逻辑）
                double indexRange_x = -1975;
                double indexRange_y = 43;
                double start_point[3] = {indexRange_x + _cached_local_x, indexRange_y + _cached_local_y, _cached_abs_alt_m};
                double end_point[3]   = {indexRange_x + target_local_x,  indexRange_y + target_local_y,  target_fixed_alt_m};
                Rect3d search_rect(target_local_x + indexRange_x, indexRange_y + target_local_y, target_fixed_alt_m,
                                  indexRange_x + _cached_local_x, indexRange_y + _cached_local_y, _cached_abs_alt_m);

                if(_rtree->Intersect3d(search_rect.min, search_rect.max, IntersectCallback)){
                    density_values[index] = 0.0f; // 0 代表被遮挡
                    continue;
                }

                // e. 坐标转换和能量计算（与原函数相同）
                float xyz[3]{};
                float thephiR[3]{};
                ct.lon2xyz(_pos.lon, _pos.lat, _pos.alt, pos_list[index].lon, pos_list[index].lat, pos_list[index].alt, xyz);
                ct.trans_axes(xyz, _base_theta, _base_phi);
                ct.xyz2tpr(xyz, thephiR);

                float minradxtheta = 0;
                float minradxphi = 0;
                float radxR = thephiR[2] / 1000.0f;

                // 计算与雷达扫描范围的关系
                float radxtheta = thephiR[0];
                float radxphi = thephiR[1];
                radxtheta = ct.normalize_angle(radxtheta);

                // 找到最小的theta
                if (radxtheta < _params.mintheta) {
                    minradxtheta = _params.mintheta - radxtheta;
                } else if ((radxtheta >= _params.mintheta) && (radxtheta <= _params.maxtheta)) {
                    minradxtheta = 0;
                } else if (radxtheta > _params.maxtheta) {
                    minradxtheta = radxtheta - _params.maxtheta;
                }

                // 找到最小的phi
                if (radxphi < _params.minphi) {
                    minradxphi = _params.minphi - radxphi;
                } else if ((radxphi >= _params.minphi) && (radxphi <= _params.maxphi)) {
                    minradxphi = 0;
                } else if (radxphi > _params.maxphi) {
                    minradxphi = radxphi - _params.maxphi;
                }

                float tpr[3];
                if (minradxtheta >= 90 || minradxphi >= 90) {
                    tpr[0] = 0;
                    tpr[1] = 0;
                    tpr[2] = -1;
                } else {
                    tpr[0] = minradxtheta;
                    tpr[1] = minradxphi;
                    tpr[2] = radxR;
                }
                density_values[index] = PointPowerDensity(tpr);
            }
        }
        proj_destroy(proj_to_wgs84);
    } // 并行区域结束

    time_t end_time = clock();
    std::cout << "--> Parallel computation finished in " << (end_time - start_time) << " ms." << std::endl;

    // 6. 颜色映射和PNG输出（与原函数相同）
    std::cout << "--> Mapping densities to colors..." << std::endl;
    float min_visible_density = 1e12, max_visible_density = 0.0f;
    size_t visible_count = 0;

    for (size_t i = 0; i < n; ++i) {
        if (density_values[i] > 0) { // 只在可见点中寻找最大最小值
            visible_count++;
            if (density_values[i] < min_visible_density) min_visible_density = density_values[i];
            if (density_values[i] > max_visible_density) max_visible_density = density_values[i];
        }
    }

    // 检查是否有足够的可见点
    if (visible_count == 0) {
        std::cerr << "Warning: No visible points found at height " << fixed_height << "m. All points are occluded." << std::endl;
    } else {
        std::cout << "--> Found " << visible_count << " visible points out of " << n << " total points ("
                  << (100.0 * visible_count / n) << "% visibility)" << std::endl;
    }

    if (min_visible_density < 1e-6f) {
        min_visible_density = 1e-6f;
    }

    // 创建PNG图像
    png_bytep* rowPointers = (png_bytep*)malloc(steps_y * sizeof(png_bytep));
    for (int i = 0; i < steps_y; i++)
        rowPointers[i] = (png_bytep)malloc(steps_x * 4 * sizeof(png_byte));

    for (int i = 0; i < steps_y; i++) {
        for (int j = 0; j < steps_x; j++) {
            size_t index = static_cast<size_t>(i) * steps_x + j;
            float density = density_values[index];
            Color color;
            color = mapDensityToColorLog_v2(density, min_visible_density, max_visible_density);
            int ii = steps_y - 1 - i; // 垂直翻转以匹配图像坐标系
            rowPointers[ii][j * 4 + 0] = color.r;
            rowPointers[ii][j * 4 + 1] = color.g;
            rowPointers[ii][j * 4 + 2] = color.b;
            rowPointers[ii][j * 4 + 3] = 255; // Alpha通道 (不透明)
        }
    }

    // 7. 保存PNG图像
    std::ostringstream filename_stream;
    filename_stream << "capable_fixed_height_power_density_" << static_cast<int>(fixed_height) << "m.png";
    std::string output_filename = filename_stream.str();

    FILE *fp = fopen(output_filename.c_str(), "wb");
    if (!fp) {
        std::cerr << "Error: Could not open file " << output_filename << " for writing" << std::endl;
        // 清理内存
        for (int i = 0; i < steps_y; i++) {
            free(rowPointers[i]);
        }
        free(rowPointers);
        delete[] density_values;
        delete[] pos_list;
        return;
    }

    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    png_infop info_ptr = png_create_info_struct(png_ptr);
    png_init_io(png_ptr, fp);
    png_set_IHDR(png_ptr, info_ptr, steps_x, steps_y, 8, PNG_COLOR_TYPE_RGB_ALPHA,
                 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    png_write_info(png_ptr, info_ptr);
    png_write_image(png_ptr, rowPointers);
    png_write_end(png_ptr, NULL);
    fclose(fp);
    png_destroy_write_struct(&png_ptr, &info_ptr);

    std::cout << "--> PNG image saved as: " << output_filename << std::endl;

    // 8. 清理内存
    for (int i = 0; i < steps_y; i++) {
        free(rowPointers[i]);
    }
    free(rowPointers);
    delete[] density_values;
    delete[] pos_list;

    std::cout << "--> Fixed-height power density calculation completed." << std::endl;
}

// 【新增函数】批量计算多个点的功率密度（Vec3d版本）
// 使用OpenMP并行，每个线程独立创建PROJ_Manager和caltools，避免重复创建开销
float Radar::CalculateSinglePointPowerDensity(const Vec3d& target_pos_decimal, caltools* ct) {
    // 将Vec3d转换为Position结构
    Position target_pos;
    DecimalToDMS(target_pos.lon, target_pos_decimal.x);
    DecimalToDMS(target_pos.lat, target_pos_decimal.y);
    target_pos.alt = target_pos_decimal.z;

    // 坐标转换和能量计算
    float xyz[3]{};
    float thephiR[3]{};

    // 使用绝对高度进行计算
    ct->lon2xyz(_pos.lon, _pos.lat, _cached_abs_alt_m,
              target_pos.lon, target_pos.lat, target_pos.alt, xyz);
    ct->trans_axes(xyz, _base_theta, _base_phi);
    ct->xyz2tpr(xyz, thephiR);

    float minradxtheta = 0;
    float minradxphi = 0;
    float radxR = thephiR[2] / 1000.0f;

    float radxtheta = thephiR[0];
    float radxphi = thephiR[1];
    radxtheta = ct->normalize_angle(radxtheta);

    // 找到theta方向上的最小角度差
    if (radxtheta < _params.mintheta) {
        minradxtheta = _params.mintheta - radxtheta;
    } else if ((radxtheta >= _params.mintheta) && (radxtheta <= _params.maxtheta)) {
        minradxtheta = 0;
    } else if (radxtheta > _params.maxtheta) {
        minradxtheta = radxtheta - _params.maxtheta;
    }

    // 找到phi方向上的最小角度差
    if (radxphi < _params.minphi) {
        minradxphi = _params.minphi - radxphi;
    } else if ((radxphi >= _params.minphi) && (radxphi <= _params.maxphi)) {
        minradxphi = 0;
    } else if (radxphi > _params.maxphi) {
        minradxphi = radxphi - _params.maxphi;
    }

    // 构造球坐标参数
    float tpr[3];
    if (minradxtheta >= 90 || minradxphi >= 90) {
        tpr[0] = 0;
        tpr[1] = 0;
        tpr[2] = -1;
    } else {
        tpr[0] = minradxtheta;
        tpr[1] = minradxphi;
        tpr[2] = radxR;
    }

    return PointPowerDensity(tpr);
}

void Radar::CalculateBatchPowerDensity(
    const std::vector<Vec3d>& target_positions,
    std::vector<float>& power_densities
) {
    size_t n = target_positions.size();
    power_densities.resize(n);

    #pragma omp parallel
    {
        // 每个线程创建一次，复用于该线程的所有点
        caltools ct;

        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < n; i++) {
            power_densities[i] = CalculateSinglePointPowerDensity(target_positions[i], &ct);
        }
    }
}


void Radar::CapablePowerDensity(int range, int level)
{
    // 1. 【新增】创建网格
    //    _pos (雷达位置) 和 _params.range (探测距离) 已在 Move() 和构造函数中设置好
    std::cout << "--> Creating grid for calculation..." << std::endl;
    auto point_list = CreateGrid(_pos, range, level); // 假设 level=16, 可后续改为参数
    auto pos_list = point_list.first;
    auto n = point_list.second;

    if (pos_list == nullptr || n == 0) {
        std::cout << "Error: CreateGrid returned a null pointer!" << std::endl;
        return;
    }
    std::cout << "--> Grid created with " << n << " points." << std::endl;

    // 3. 分配用于存储最终结果的内存
    float *capvalues = new float[n];
    std::cout << "--> Starting parallel computation..." << std::endl;
    time_t start_time = clock();
    
    // 4. 调用上面那个并行的、底层的计算函数
    CapablePowerDensity(capvalues, pos_list, n);

    time_t end_time = clock();
    std::cout << "--> Parallel computation finished in " << (end_time - start_time) << " ms." << std::endl;

    // 5. 将计算结果批量写入数据库
    std::cout << "--> Writing " << n << " results to database..." << std::endl;
    _writer->Write(pos_list, capvalues, n, _timestamp);
    std::cout << "--> Write operation completed." << std::endl;
    delete[] capvalues;
    delete[] pos_list;
}

// /**
//  * 它计算的不是雷达某一瞬间的功率密度，而是当雷达天线在其允许的扫描范围内扫视时，
//  * 每个网格点可能接收到的最大功率密度。它生成的是一张“最优探测能力图”。
//  */
// void Radar::CapablePowerDensity(float out[], const Position *pos_list, size_t n, float radius[])
// {
//     int p_n, p_id;
//     MPI_Comm_size(MPI_COMM_WORLD, &p_n);
//     MPI_Comm_rank(MPI_COMM_WORLD, &p_id);
//     auto chunk = n / p_n;
//     auto start = p_id * chunk;
//     auto end = p_id < p_n - 1 ? (p_id + 1) * chunk : n;
//     auto local_n = end - start;
//     caltools ct;
//     float(*tpr_local)[3] = (float(*)[3])malloc(sizeof(float) * local_n * 3);
//     time_t start_time = clock();
//     for (int i = 0; i < local_n; i++)
//     {
//         int glocal_i = i + start;
//         const auto& target_pos = pos_list[glocal_i];
//         if (_rtree != nullptr && isOccluded(target_pos))
//         {
//             tpr_local[i][0] = 0;
//             tpr_local[i][1] = 0;
//             tpr_local[i][2] = -1;
//             continue;
//         }
//         float xyz[3]{};
//         float thephiR[3]{};
//         ct.lon2xyz(_pos.lon, _pos.lat, _pos.alt, pos_list[glocal_i].lon, pos_list[glocal_i].lat, pos_list[glocal_i].alt, xyz);
//         ct.trans_axes(xyz, _base_theta, _base_phi);
//         ct.xyz2tpr(xyz, thephiR);
//         //现在thephiR是p点与载具初始方向的夹角
//         float minradxtheta = 0;
//         float minradxphi = 0;
//         float radxR = thephiR[2] / 1000; //米转换为千米
//         if (radius != NULL)
//         {
//             radius[glocal_i] = radxR;
//         }
//         //先计算每点与载具角度的差值，再看其是否在旋转范围内
//         float radxtheta = thephiR[0];
//         float radxphi = thephiR[1];
//         radxtheta = ct.normalize_angle(radxtheta);
//         //找到最小的theta
//         if (radxtheta < _params.mintheta)
//         {
//             minradxtheta = _params.mintheta - radxtheta;
//         }
//         else if ((radxtheta >= _params.mintheta) && (radxtheta <= _params.maxtheta))
//         {
//             minradxtheta = 0;
//         }
//         else if (radxtheta > _params.maxtheta)
//         {
//             minradxtheta = radxtheta - _params.maxtheta;
//         }
//         //找到最小的phi
//         if (radxphi < _params.minphi)
//         {
//             minradxphi = _params.minphi - radxphi;
//         }
//         else if ((radxphi >= _params.minphi) && (radxphi <= _params.maxphi))
//         {
//             minradxphi = 0;
//         }
//         else if (radxphi > _params.maxphi)
//         {
//             minradxphi = radxphi - _params.maxphi;

//         }

//         if (minradxtheta >= 90 || minradxphi >= 90)
//         {
//             tpr_local[i][0] = 0;
//             tpr_local[i][1] = 0;
//             tpr_local[i][2] = -1;
//         }
//         else
//         {
//             tpr_local[i][0] = minradxtheta;
//             tpr_local[i][1] = minradxphi;
//             tpr_local[i][2] = radxR;
//         }
//     }
//     //以上是数据转换，以下是模型计算
//     for (int i = 0; i < local_n; i++)
//     {
//         out[start + i] = PointPowerDensity(tpr_local[i]);
//     }
//     time_t end_time = clock();
// #ifdef DEBUG
//     if (p_id == 0)
//     {
//         cout << "网格数：" << n << "，计算用时：";
//         cout << (end_time - start_time) << endl;
//     }
// #endif
//     free(tpr_local);
// }
// void Radar::CapablePowerDensity(const Position *pos_list, size_t n)
// {
//     // TODO: refactor to remove duplicates
//     int p_n, p_id;
//     MPI_Comm_size(MPI_COMM_WORLD, &p_n);
//     MPI_Comm_rank(MPI_COMM_WORLD, &p_id);
//     auto chunk = n / p_n;
//     auto start = p_id * chunk;
//     auto end = p_id < p_n - 1 ? (p_id + 1) * chunk : n;
//     auto local_n = end - start;

//     float *capvalues = (float *)malloc(sizeof(float) * n);
//     CapablePowerDensity(capvalues, pos_list, n);
//     _writer->Write(pos_list + start, capvalues + start, local_n, _timestamp);
//     free(capvalues);
// }

void Radar::CapableRegion(const Position *pos_list, size_t n, float rcs)
{
    // int p_n, p_id;
    // MPI_Comm_size(MPI_COMM_WORLD, &p_n);
    // MPI_Comm_rank(MPI_COMM_WORLD, &p_id);
    // auto chunk = n / p_n;
    // auto start = p_id * chunk;
    // auto end = p_id < p_n - 1 ? (p_id + 1) * chunk : n;
    // auto local_n = end - start;

    // float *capvalues = (float *)malloc(sizeof(float) * n);
    // float *radius = (float *)malloc(sizeof(float) * n);

    // CapablePowerDensity(capvalues, pos_list, n, radius);
    // auto range_km = _params.range / 1000;
    // float ref_tpr[3] = {0, 0, range_km};
    // // k = density / r^2
    // float ref_k = PointPowerDensity(ref_tpr) * rcs / 7.5 / (range_km * range_km);
    // for (size_t i = start; i < end; i++)
    // {
    //     capvalues[i] = (capvalues[i] / (radius[i] * radius[i]) < ref_k) ? 0 : 1;
    // }
    // _writer->Write(pos_list + start, capvalues + start, local_n, _timestamp);
    // free(capvalues);
    // free(radius);
}
// 添加一个辅助函数，用于将度分秒转换为十进制度

// TODO 之后需要进行改善
// double haversineDistance(float lon1, float lat1, float lon2, float lat2) {
//     double dLat = (lat2 - lat1) * PI / 180.0;
//     double dLon = (lon2 - lon1) * PI / 180.0;

//     lat1 = lat1 * PI / 180.0;
//     lat2 = lat2 * PI / 180.0;

//     double a = pow(sin(dLat / 2), 2) + pow(sin(dLon / 2), 2) * cos(lat1) * cos(lat2);
//     double c = 2 * asin(sqrt(a));
//     return EARTH_RADIUS_M * c;
// }
// // 低精度、高速度的距离估算函数
// double fastApproximateDistance(float lon1, float lat1, float lon2, float lat2) {
//     // 将纬度从度转换为弧度
//     double lat1_rad = lat1 * PI / 180.0;
    
//     // 计算经纬度差（度）
//     double dLon = lon2 - lon1;
//     double dLat = lat2 - lat1;

//     // 将经纬度差从度转换为米（近似值）
//     // 1度纬度约等于 111132.954 米
//     // 1度经度约等于 111319.488 * cos(latitude) 米
//     double dy = dLat * 111132.954;
//     double dx = dLon * 111319.488 * cos(lat1_rad);

//     // 使用勾股定理
//     return sqrt(dx * dx + dy * dy);
// }

// 5. 实现新的 isOccluded 核心逻辑函数
bool Radar::isOccluded(const Position& target_pos,PROJ_Manager* proj_manager) const {
    if(_rtree == nullptr || _cached_abs_alt_m < -9000.0 || proj_manager == nullptr){
        return false;
    }
        // 1. 【保留】准备终点（目标）的坐标
    double target_lon_deg = DMSToDecimal(target_pos.lon);
    double target_lat_deg = DMSToDecimal(target_pos.lat);
    double target_abs_alt_m = target_pos.alt;

    // 2. 【保留】只对终点执行坐标转换
    double min_x = 835000;
    double min_y = 815500;
    double indexRange_x = -1975;
    double indexRange_y = 43;

    PJ_COORD coord_end = proj_coord(target_lat_deg, target_lon_deg, 0, 0);
    PJ_COORD result_end = proj_trans(proj_manager->proj_from, PJ_FWD, coord_end);

    double target_local_x = result_end.xy.y - min_x;
    double target_local_y = result_end.xy.x - min_y;
    // 3. 【修改】准备查询坐标，起点直接使用缓存好的值
    double start_point[3] = {indexRange_x + _cached_local_x, indexRange_y + _cached_local_y, _cached_abs_alt_m};
    double end_point[3]   = {indexRange_x + target_local_x,  indexRange_y + target_local_y,  target_abs_alt_m};
      
    // 4. 执行查询 (不变)
    return _rtree->Intersect3d(start_point, end_point, IntersectCallback);

}
