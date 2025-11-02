#pragma once
#include "hiradar/radar.hpp"

inline Radar *ANAPG66(Position &pos, float base_phi, float base_theta)
{
    RadarParams params = {"AN/APG-66", RadarType::SEW, 1, 0.48, 10000, 32, 20, -30, 30, -30, 30, 46000};
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