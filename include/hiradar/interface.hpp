#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <cstring>
#include <map>

#define PRECISION 1000 //角度精度
// #define PI 3.1415926
#define EARTH_RADIUS 6371.193 //地球半径 KM

//以度秒分标识经纬度，以米标识高程
struct Position
{
    float lon[3] = {0};
    float lat[3] = {0};
    float alt = 0;
};

class IPowerDensityWriter
{
public:
    virtual void Write(const Position *pos_list, float *values, size_t n, float timestamp) = 0;
    virtual ~IPowerDensityWriter();
};

inline void PrintDegree(std::ostream &out, float degree[3])
{
    out << degree[0] << "°" << degree[1] << "'" << degree[2] << "\"";
}

class StdoutDensityWriter : public IPowerDensityWriter
{
public:
    virtual void Write(const Position *pos_list, float *values, size_t n, float timestamp)
    {
        for (size_t i = 0; i < n; i++)
        {
            auto pos = pos_list[i];
            auto val = values[i];
            std::cout << "[" << timestamp << "],";
            PrintDegree(std::cout, pos.lon);
            std::cout << ",";
            PrintDegree(std::cout, pos.lat);
            std::cout << "," << pos.alt << " " << val << std::endl;
        }
    }
};

//这是一个抽象类，用于表示雷达，雷达可以移动，旋转，计算功率密度，计算可覆盖区域

class IRadar
{
public:
    virtual ~IRadar();

    virtual void Move(const Position &pos, float base_phi, float base_theta, float timestamp) {} // 移动
    virtual void Rotate(float phi, float theta, float timestamp) = 0; // 旋转
    virtual void PowerDensity(const Position *pos_list, size_t n) = 0; // 计算n个点的能量密度 
    virtual void CapablePowerDensity(int range, int level) = 0; // n个点的最大能量密度
    virtual void CapableRegion(const Position *pos_list, size_t n, float rcs) = 0; // 指定rcs所能探测到的网格

    virtual void BindWriter(IPowerDensityWriter *writer) { _writer = writer; }; //绑定写入器

    // getters
    IPowerDensityWriter *GetWriter() const { return _writer; }
    virtual Position GetPosition() const = 0;

protected:
    IPowerDensityWriter *_writer = NULL;
};

//这是啥？
const float RCSMiddle = 0;
const float RCSSmall = 0;
const float RCSLarge = 0;

enum RadarType
{
    SEW,    //缝阵
    PHASED, //相控阵
    YAGI    //八木
};
struct RadarParams
{
    std::string name;              //雷达名称
    RadarType type;                //雷达类型
    int N;                         //阵元个数
    float d;                       //阵列间距或y方向尺寸
    float f;                       //工作频率，单位为Mhz
    float G;                       //天线增益。单位为dB
    float Pt;                      //峰值功率,单位为kW
    float mintheta;                //最小的theta，逆时针为正方向
    float maxtheta;                //最大的theta
    float minphi;                  //最小的phi
    float maxphi;                  //最大的phi
    float range;                   //RCS=7.5下雷达作用距离(米为单位）
};
