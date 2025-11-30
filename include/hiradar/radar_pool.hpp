#pragma once
#include "hiradar/radar.hpp"

inline Radar *ANAPG66(Position &pos, float base_phi, float base_theta)
{
    // 1. 定义一个RadarParams结构体，并用硬编码的数值来初始化它
    //    这些数值就是AN/APG-66这款雷达独有的、固定的物理特性
    RadarParams params = {"AN/APG-66", RadarType::SEW, 1, 0.48, 10000, 32, 20, -30, 30, 0, 30, 46000};
    // 2. 使用 new 关键字创建一个通用的 Radar 对象实例
    //    在创建时，将上面定义好的、包含特定型号参数的 params 结构体，
    //    以及初始的位置和姿态信息，传递给 Radar 类的构造函数
    return new Radar(params, pos, base_phi, base_theta);
}

inline Radar *ANAPS115(Position &pos, float base_phi, float base_theta)
{
    RadarParams params = {"AN/APS115", RadarType::SEW, 1, 0.66, 10000, 35, 143, -45, 45, -20, 10, 200000};
    return new Radar(params, pos, base_phi, base_theta);
}

inline Radar *ANAPY12(Position &pos, float base_phi, float base_theta)
{
    RadarParams params = {"ANAPY12", RadarType::SEW, 1, 1.3, 3000, 29, 500, -180, 180, -30, 30, 445000};
    return new Radar(params, pos, base_phi, base_theta);
}

inline Radar *ANSPS49(Position &pos, float base_phi, float base_theta)
{
    RadarParams params = {"ANSPS49", RadarType::SEW, 1, 4.3, 1000, 28.5, 360, -180, 180, -23.5, 23.5, 463000};
    return new Radar(params, pos, base_phi, base_theta);
}

inline Radar *SPY1D(Position &pos, float base_phi, float base_theta)
{
    RadarParams params = {"SPY1D", RadarType::PHASED, 4480, 3.5, 3000, 42, 5000, -180, 180, 0, 90, 324100};
    return new Radar(params, pos, base_phi, base_theta);
}

inline Radar *ANMPQ65(Position &pos, float base_phi, float base_theta)
{
    RadarParams params = {"ANMPQ65", RadarType::PHASED, 5161, 2.44, 5000, 39, 200, -60, 60, 22.5, 85, 170000};
    return new Radar(params, pos, base_phi, base_theta);
}

inline Radar *ANAPS145(Position &pos, float base_phi, float base_theta)
{
    RadarParams params = {"ANAPS145", RadarType::YAGI, 1, 1, 450, 32, 1000, -180, 180, 0, 0, 270000};
    return new Radar(params, pos, base_phi, base_theta);
}
inline Radar *CityGuardRadar(Position &pos, float base_phi, float base_theta){
    RadarParams params = {
        "CityGuardRadar",    // 1. 雷达名称
        RadarType::PHASED,   // 2. 雷达类型 (相控阵雷达适合城市监控)
        256,                 // 3. 阵元个数 (可以少一些)
        0.006,                // 4. 阵元间距 (米，可以小一些以适应高频率)
        24000,               // 5. 工作频率 (MHz, K波段频率高，适合近距离高分辨率)
        25,                  // 6. 天线增益 (dB)
        0.001,                  // 7. 峰值功率 (kW，近距离不需要太高功率)
        -180, 180,             // 8. 方位角扫描范围 (提供一个较宽的监控扇区)
        -90, 0,             // 9. 俯仰角扫描范围 (适合从高处向下监控)
        300                  // 10.【关键】作用距离 (米)，设置为800米，完美覆盖500x500的区域
    };
    return new Radar(params, pos, base_phi, base_theta);
}