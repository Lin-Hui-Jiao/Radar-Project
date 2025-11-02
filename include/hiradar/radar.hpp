#pragma once
#include <assert.h>
#include "hiradar/interface.hpp"
#include "hiradar/caltools.h"
// #include "hiradar/TerrainGrid.hpp"
#include "hiradar/query3dRtree.h"
#include "hiradar/RTree.h"
#include "hiradar/occlusion_utils.h"


// general radar for both sew and phased
class Radar : public IRadar
{
public:
    //构造函数，初始化雷达参数，雷达位置，雷达自身俯仰角，雷达自身方位角
    Radar(RadarParams params, Position &pos, float base_phi = 0, float base_theta = 0)
        : _params(params), _etheta_table(nullptr), _ephi_table(nullptr), _rtree(nullptr),
          _pos(pos), _base_phi(base_phi), _base_theta(base_theta), _phi(0), _theta(0), _timestamp(0),
          _cached_local_x(0.0), _cached_local_y(0.0), _cached_abs_alt_m(-9999.0)
    {
        // 关键优化：为方向性系数创建预计算查找表 (Pre-computed Lookup Table)
        // 这一步对应了论文中提到的“预缓存加速策略”

        // 为方位角(theta)的方向性系数分配内存，大小为 360 * 1000，1000是精度
        _etheta_table = (float *)malloc(sizeof(float) * 360 * PRECISION);
        // 为方位角(phi)的方向性系数分配内存，大小为 360 * 1000，1000是精度
        _ephi_table = (float *)malloc(sizeof(float) * 360 * PRECISION);

        // 调用 caltools 中的函数，预先计算并填充方位角查找表
        caltools().make_angle_table(_params.type, _params.N, _params.d, _params.f, _etheta_table,true);//true是计算ethetatable
        // 调用 caltools 中的函数，预先计算并填充方位角查找表 
        caltools().make_angle_table(_params.type, _params.N, _params.d, _params.f, _ephi_table,false);//false是计算ephitable
    }

    ~Radar() {
        free(_etheta_table);
        free(_ephi_table);
        }

    // 实现IRadar接口中的Move函数，用于更新雷达载具的位置和基础姿态
    virtual void Move(const Position &pos, float base_phi, float base_theta, float timestamp)
    {
        _pos = pos;
        _base_phi = base_phi;
        _base_theta = base_theta;
        _timestamp = timestamp;

        if (_rtree == nullptr) {
            // 如果没有绑定，可以将缓存设为无效值或打印警告
            _cached_abs_alt_m = -9999.0; // 例如，设为一个无效高度
            return;
        }
         // 将雷达位置（固定点）的转换逻辑从 isOccluded 移动到这里
        double radar_lon_deg = DMSToDecimal(_pos.lon);
        double radar_lat_deg = DMSToDecimal(_pos.lat);
        double radar_alt_m_relative = _pos.alt; //注意雷达的高度是相对于地面的高度，在开始代码执行的时候，传入的参数一定要注意

        // 坐标系常量
        double min_x = 835000;
        double min_y = 815500;
        double indexRange_x = -1975;
        double indexRange_y = 43;

        // 执行PROJ转换
        PJ_CONTEXT *context = proj_context_create();
        PJ *proj_from = proj_create_crs_to_crs(context, "EPSG:4326", "EPSG:2326", NULL);
        PJ_COORD coord_start = proj_coord(radar_lat_deg, radar_lon_deg, 0, 0);
        PJ_COORD result_start = proj_trans(proj_from, PJ_FWD, coord_start);
        proj_destroy(proj_from);
        proj_context_destroy(context);

        // 计算局部坐标并存入缓存
        _cached_local_x = result_start.xy.y - min_x;
        _cached_local_y = result_start.xy.x - min_y;

        // 查询地表高度并计算绝对高度，然后存入缓存
        double terrain_alt_m = _rtree->Getheight3d(indexRange_x + _cached_local_x, indexRange_y + _cached_local_y, HeightCallback);
        _cached_abs_alt_m = terrain_alt_m + radar_alt_m_relative;
        cout << "雷达的绝对海拔是: " << _cached_abs_alt_m << endl;
        cout << "雷达的相对海拔是: " << radar_alt_m_relative << endl;
    }

    // 实现IRadar接口中的Rotate函数，用于更新雷达自身姿态
    virtual void Rotate(float phi, float theta, float timestamp)
    {
        _phi = phi;
        _theta = theta;
        _timestamp = timestamp;
    }

    // 声明核心计算函数，它们的具体实现代码在对应的 radar.cpp 文件中
    
    // 计算给定点列表在当前雷达姿态下的功率密度
    virtual void PowerDensity(const Position *pos_list, size_t n);

    // 计算给定点列表在雷达整个扫描范围内的最大可能功率密度
    // virtual void CapablePowerDensity(const Position *pos_list, size_t n);

    virtual void CapablePowerDensity(int range, int level);

    // 计算给定点列表中，哪些区域对于一个特定RCS值的目标是可探测的
    virtual void CapableRegion(const Position *pos_list, size_t n, float rcs);

    void CapableGroundPowerDensity(int range,double step);

    void CapableFixedHeightPowerDensity(int range, double step, double fixed_height);

    // 实现IRadar接口中的GetPosition函数，用于获取雷达当前位置
    virtual Position GetPosition() const { return _pos; }

    void BindRtree(RTree3d* rtree) { _rtree = rtree; }

    // 【新增】批量计算多个点的功率密度（Vec3d版本）
    // 输入：target_positions - 目标点列表（Vec3d格式：x=经度(十进制度), y=纬度(十进制度), z=绝对高度(米)）
    // 输出：power_densities - 对应的功率密度（W/m²），如果被遮挡或超出扫描范围则为0
    // 注意：
    // 1. 内部使用OpenMP并行，每个线程独立创建PROJ_Manager和caltools
    // 2. 不进行遮挡检测，调用者需要预先筛选可见点
    // 3. power_densities会被自动resize为与target_positions相同大小
    void CalculateBatchPowerDensity(
        const std::vector<Vec3d>& target_positions,
        std::vector<float>& power_densities
    );

    // 【新增】计算单个点的功率密度（Vec3d版本，供OpenMP并行使用）
    // 输入：target_pos - 目标点（Vec3d格式：x=经度(十进制度), y=纬度(十进制度), z=绝对高度(米)）
    //      ct - 线程局部的caltools对象指针（由调用者提供，避免重复创建）
    // 输出：功率密度（W/m²），如果超出扫描范围则为0
    // 注意：不进行遮挡检测，调用者需要预先判断可见性
    float CalculateSinglePointPowerDensity(const Vec3d& target_pos, caltools* ct);

    double _cached_abs_alt_m;

private:
    RadarParams _params;
    float *_etheta_table;
    float *_ephi_table;
    RTree3d * _rtree = nullptr;

    // dynamics
    Position _pos;     //雷达位置
    float _base_phi;   //载具自身俯仰角
    float _base_theta; //载具自身方位角
    float _phi;        //雷达自身俯仰角
    float _theta;      //雷达自身方位角
    float _timestamp;  // unix time, can be relative
    double _cached_local_x;
    double _cached_local_y;

    // 计算单个点的功率密度，是内部使用的辅助函数
    float PointPowerDensity(float tpr[3]);
    void CapablePowerDensity(float out[], const Position *pos_list, size_t n, float radius[] = NULL);
    
    bool isOccluded(const Position& target_pos,PROJ_Manager* proj_manager) const;
};
