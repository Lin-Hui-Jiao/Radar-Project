#include "hiradar/radar.hpp"
#include <cmath>
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
void Radar::PowerDensity(const Position *pos_list, size_t n)
{
    int p_n, p_id;
    MPI_Comm_size(MPI_COMM_WORLD, &p_n);
    MPI_Comm_rank(MPI_COMM_WORLD, &p_id);
    auto chunk = n / p_n;
    auto start = p_id * chunk;
    auto end = p_id < p_n - 1 ? (p_id + 1) * chunk : n;
    auto local_n = end - start;
    float *values = (float *)malloc(sizeof(float) * n);
    caltools ct;
    float(*tpr)[3] = (float(*)[3])malloc(sizeof(float) * n * 3);
    time_t start_time = clock();
    //坐标变换流水线
    for (int i = start; i < end; i++)
    {
        float xyz[3]{};
        float thephiR[3]{};
        //全局坐标系转为雷达所在位置为原点，正北方向为x轴
        ct.lon2xyz(_pos.lon, _pos.lat, _pos.alt, pos_list[i].lon, pos_list[i].lat, pos_list[i].alt, xyz);
        //变换坐标系之后，点在新坐标系中的位置(x,y,z)
        ct.trans_axes(xyz, _base_theta + _theta, _base_phi + _phi);
        // 笛卡尔坐标系转球坐标系
        ct.xyz2tpr(xyz, thephiR);
        //现在thephiR是p点与载具加上雷达方向的夹角
        float radxR = thephiR[2] / 1000; //米转换为千米

        tpr[i][0] = thephiR[0];
        tpr[i][1] = thephiR[1];
        tpr[i][2] = radxR;
    }

    //这里是错误
    for (int i = 0; i < n; i++)
    {
        values[i] = PointPowerDensity(tpr[i]);
    }

    time_t end_time = clock();
#ifdef DEBUG
    if (p_id == 0)
    {
        cout << "网格数：" << n << "，计算用时：";
        cout << (end_time - start_time) << endl;
    }
#endif

    free(tpr);
    _writer->Write(pos_list + start, values + start, local_n, _timestamp);//并发的写入到数据库中
    free(values);
}
/**
 * 它计算的不是雷达某一瞬间的功率密度，而是当雷达天线在其允许的扫描范围内扫视时，
 * 每个网格点可能接收到的最大功率密度。它生成的是一张“最优探测能力图”。
 */
void Radar::CapablePowerDensity(float out[], const Position *pos_list, size_t n, float radius[])
{
    int p_n, p_id;
    MPI_Comm_size(MPI_COMM_WORLD, &p_n);
    MPI_Comm_rank(MPI_COMM_WORLD, &p_id);
    auto chunk = n / p_n;
    auto start = p_id * chunk;
    auto end = p_id < p_n - 1 ? (p_id + 1) * chunk : n;
    auto local_n = end - start;
    caltools ct;
    float(*tpr_local)[3] = (float(*)[3])malloc(sizeof(float) * local_n * 3);
    time_t start_time = clock();
    for (int i = 0; i < local_n; i++)
    {
        int glocal_i = i + start;
        const auto& target_pos = pos_list[glocal_i];
        if (_terrain_grid != nullptr && isOccluded(target_pos))
        {
            tpr_local[i][0] = 0;
            tpr_local[i][1] = 0;
            tpr_local[i][2] = -1;
            continue;
        }
        float xyz[3]{};
        float thephiR[3]{};
        ct.lon2xyz(_pos.lon, _pos.lat, _pos.alt, pos_list[glocal_i].lon, pos_list[glocal_i].lat, pos_list[glocal_i].alt, xyz);
        ct.trans_axes(xyz, _base_theta, _base_phi);
        ct.xyz2tpr(xyz, thephiR);
        //现在thephiR是p点与载具初始方向的夹角
        float minradxtheta = 0;
        float minradxphi = 0;
        float radxR = thephiR[2] / 1000; //米转换为千米
        if (radius != NULL)
        {
            radius[glocal_i] = radxR;
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

        if (minradxtheta >= 90 || minradxphi >= 90)
        {
            tpr_local[i][0] = 0;
            tpr_local[i][1] = 0;
            tpr_local[i][2] = -1;
        }
        else
        {
            tpr_local[i][0] = minradxtheta;
            tpr_local[i][1] = minradxphi;
            tpr_local[i][2] = radxR;
        }
    }
    //以上是数据转换，以下是模型计算
    for (int i = 0; i < local_n; i++)
    {
        out[start + i] = PointPowerDensity(tpr_local[i]);
    }
    time_t end_time = clock();
#ifdef DEBUG
    if (p_id == 0)
    {
        cout << "网格数：" << n << "，计算用时：";
        cout << (end_time - start_time) << endl;
    }
#endif
    free(tpr_local);
}
void Radar::CapablePowerDensity(const Position *pos_list, size_t n)
{
    // TODO: refactor to remove duplicates
    int p_n, p_id;
    MPI_Comm_size(MPI_COMM_WORLD, &p_n);
    MPI_Comm_rank(MPI_COMM_WORLD, &p_id);
    auto chunk = n / p_n;
    auto start = p_id * chunk;
    auto end = p_id < p_n - 1 ? (p_id + 1) * chunk : n;
    auto local_n = end - start;

    float *capvalues = (float *)malloc(sizeof(float) * n);
    CapablePowerDensity(capvalues, pos_list, n);
    _writer->Write(pos_list + start, capvalues + start, local_n, _timestamp);
    free(capvalues);
}

void Radar::CapableRegion(const Position *pos_list, size_t n, float rcs)
{
    int p_n, p_id;
    MPI_Comm_size(MPI_COMM_WORLD, &p_n);
    MPI_Comm_rank(MPI_COMM_WORLD, &p_id);
    auto chunk = n / p_n;
    auto start = p_id * chunk;
    auto end = p_id < p_n - 1 ? (p_id + 1) * chunk : n;
    auto local_n = end - start;

    float *capvalues = (float *)malloc(sizeof(float) * n);
    float *radius = (float *)malloc(sizeof(float) * n);

    CapablePowerDensity(capvalues, pos_list, n, radius);
    auto range_km = _params.range / 1000;
    float ref_tpr[3] = {0, 0, range_km};
    // k = density / r^2
    float ref_k = PointPowerDensity(ref_tpr) * rcs / 7.5 / (range_km * range_km);
    for (size_t i = start; i < end; i++)
    {
        capvalues[i] = (capvalues[i] / (radius[i] * radius[i]) < ref_k) ? 0 : 1;
    }
    _writer->Write(pos_list + start, capvalues + start, local_n, _timestamp);
    free(capvalues);
    free(radius);
}
// 添加一个辅助函数，用于将度分秒转换为十进制度
float dmsToDecimal(const float dms[3]) {
    return dms[0] + dms[1] / 60.0f + dms[2] / 3600.0f;
}
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
// 低精度、高速度的距离估算函数
double fastApproximateDistance(float lon1, float lat1, float lon2, float lat2) {
    // 将纬度从度转换为弧度
    double lat1_rad = lat1 * PI / 180.0;
    
    // 计算经纬度差（度）
    double dLon = lon2 - lon1;
    double dLat = lat2 - lat1;

    // 将经纬度差从度转换为米（近似值）
    // 1度纬度约等于 111132.954 米
    // 1度经度约等于 111319.488 * cos(latitude) 米
    double dy = dLat * 111132.954;
    double dx = dLon * 111319.488 * cos(lat1_rad);

    // 使用勾股定理
    return sqrt(dx * dx + dy * dy);
}

// 5. 实现新的 isOccluded 核心逻辑函数
bool Radar::isOccluded(const Position& target_pos) const {
    float radar_lon_deg = dmsToDecimal(_pos.lon);
    float radar_lat_deg = dmsToDecimal(_pos.lat);
    float radar_alt_m = _pos.alt;

    float target_lon_deg = dmsToDecimal(target_pos.lon);
    float target_lat_deg = dmsToDecimal(target_pos.lat);
    float target_alt_m = target_pos.alt;

    float distance = fastApproximateDistance(radar_lon_deg, radar_lat_deg, target_lon_deg, target_lat_deg);
    float ground_resolution = _terrain_grid->getGroundResolution();
    if (ground_resolution <= 0) {
        // 处理无效分辨率
        return false; // 或其他合适的处理
    }
    int steps = static_cast<int>(ceil(distance / ground_resolution)) + 1;

    float lon_step = (target_lon_deg - radar_lon_deg) / steps;
    float lat_step = (target_lat_deg - radar_lat_deg) / steps;
    float alt_step = (target_alt_m - radar_alt_m) / steps;

    // 从第1步开始检查（不检查雷达脚下），直到目标点前一步
    for (int i = 1; i < steps; ++i) {
        float current_lon = radar_lon_deg + i * lon_step;
        float current_lat = radar_lat_deg + i * lat_step;
        float los_alt_m = radar_alt_m + i * alt_step;

        float terrain_alt_m = _terrain_grid->getNearestElevation(current_lon, current_lat);

        // 如果地形有效(-9999.0f是无效值)且视线高度低于或等于地形高度，则被遮挡
        if (terrain_alt_m > -9000.0f && los_alt_m < terrain_alt_m) {
            return true; // 被遮挡
        }
    }
    return false; // 未被遮挡
}
