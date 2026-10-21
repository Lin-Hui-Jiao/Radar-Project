// radar_kernels.hpp — amplitude kernels for field_export, copied verbatim from the PSHC prototype
// (pshc/prototypes/rise_prototype/src/common.hpp: struct Meta, ReplayKernel, CleanKernel).
//   ReplayKernel: line-by-line port of the float32 chain of Radar::CalculateSinglePointPowerDensity
//                 (src/radar.cpp, src/caltools.cpp); bit-identical to the existing program  -> data version V-float
//   CleanKernel : the same formulas in double precision                                     -> data version V-double
// Do not edit the kernels: bit-exactness against build/bin/fixed_height_result_*.csv depends on every operation.
#pragma once
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>

// ---------------------------------------------------------------------------
// Grid + radar metadata (plain "key value" text file)
// ---------------------------------------------------------------------------
struct Meta {
    int nx = 0, ny = 0, nz = 0;            // lon, lat, alt counts
    double lon0 = 0, lat0 = 0, alt0 = 0;   // first node
    double dlon = 0, dlat = 0, dalt = 0;   // node spacing (deg, deg, m)
    double rlon = 0, rlat = 0, ralt = 0;   // radar position (deg, deg, abs alt m)
    double mlon = 0, mlat = 0;             // metres per degree used by the occlusion geometry
    // radar parameters (CityGuardRadar defaults, include/hiradar/radar_pool.hpp)
    int N = 256; double d = 0.006, f = 24000, G = 25, Pt = 0.001;
    double mintheta = -180, maxtheta = 180, minphi = -90, maxphi = 0;
    size_t count() const { return (size_t)nx * ny * nz; }
    double lon(int i) const { return lon0 + i * dlon; }
    double lat(int j) const { return lat0 + j * dlat; }
    double alt(int k) const { return alt0 + k * dalt; }
    // local metric coordinates relative to the radar (same affine frame the occlusion test uses)
    inline void local(int i, int j, int k, double& e, double& n, double& u) const {
        e = (lon(i) - rlon) * mlon; n = (lat(j) - rlat) * mlat; u = alt(k) - ralt;
    }
    inline void localOf(double lo, double la, double al, double& e, double& n, double& u) const {
        e = (lo - rlon) * mlon; n = (la - rlat) * mlat; u = al - ralt;
    }
};

// ---------------------------------------------------------------------------
// Amplitude kernel A(x): unobstructed direct-path power density.
//
// ReplayKernel is a line-by-line port of the float32 pipeline in
//   src/radar.cpp  Radar::CalculateSinglePointPowerDensity / PointPowerDensity
//   src/caltools.cpp lon2xyz / cal_distance / xyz2tpr / make_angle_table / cal_Qt
// (trans_axes is the identity for base_theta = base_phi = 0 and is omitted).
// It reproduces the float32 rounding behaviour of the generator, including the
// loss of horizontal precision caused by holding decimal degrees in float.
//
// CleanKernel evaluates the same formulas in double precision.
// ---------------------------------------------------------------------------
struct ReplayKernel {
    static constexpr double PIc = 3.1415926;      // #define PI in caltools.cpp
    static constexpr double EARTH_RADIUS = 6371.193;
    static constexpr int PRECISION = 1000;
    std::vector<float> tab;                         // 360*PRECISION entries
    float Rlon[3], Rlat[3]; float Ralt;
    int N; float d, f, G, Pt, mintheta, maxtheta, minphi, maxphi;

    static void DecimalToDMS(float dms[3], float de) {
        int D = static_cast<int>(de);
        int M = static_cast<int>((de - D) * 60);
        float S = (de - D) * 3600.0f - M * 60.0f;
        dms[0] = D; dms[1] = M; dms[2] = S;
    }
    static float rad(float dd) { return dd * PIc / 180.0; }
    static float cal_distance(float lat1, float lng1, float lat2, float lng2) {
        float a, b;
        float radLat1 = rad(lat1), radLat2 = rad(lat2);
        a = radLat1 - radLat2;
        b = rad(lng1) - rad(lng2);
        float s = 2 * std::asin(std::sqrt(std::pow(std::sin(a / 2), 2) + std::cos(radLat1) * std::cos(radLat2) * std::pow(std::sin(b / 2), 2)));
        s = s * EARTH_RADIUS;
        s = s * 1000;
        return s;
    }
    static float cal_phased_Eangle(double angle, float n, float d, float f) {
        float lamda = static_cast<float>(300) / f;
        float aaa = std::sin(angle / 180 * PIc);
        float bbb = std::sin(n * (PIc * d / lamda) * aaa);
        float ccc = n * std::sin((PIc * d / lamda) * aaa);
        return bbb / ccc;
    }
    void init(const Meta& m) {
        N = m.N; d = (float)m.d; f = (float)m.f; G = (float)m.G; Pt = (float)m.Pt;
        mintheta = (float)m.mintheta; maxtheta = (float)m.maxtheta; minphi = (float)m.minphi; maxphi = (float)m.maxphi;
        DecimalToDMS(Rlon, (float)m.rlon); DecimalToDMS(Rlat, (float)m.rlat); Ralt = (float)m.ralt;
        tab.assign((size_t)360 * PRECISION, 0.f);
        for (int i = 1; i < 360 * PRECISION; ++i) tab[i] = cal_phased_Eangle(i / static_cast<double>(PRECISION), (float)N, d, f);
        tab[0] = 1; tab[180 * PRECISION] = 0;
    }
    static float normalize_angle(float angle) {
        auto a360 = std::fmod((double)angle, 360.0);
        return a360 > 180 ? 360 - angle : angle;
    }
    // (lon, lat) decimal degrees, alt metres absolute -> W/m^2 (float32)
    float operator()(double lon, double lat, double alt) const {
        float xlon[3], xlat[3]; DecimalToDMS(xlon, (float)lon); DecimalToDMS(xlat, (float)lat);
        float xalt = (float)alt;
        // lon2xyz
        float Rlondu = Rlon[0] + Rlon[1] / 60 + Rlon[2] / 3600;
        float Rlatdu = Rlat[0] + Rlat[1] / 60 + Rlat[2] / 3600;
        float xlondu = xlon[0] + xlon[1] / 60 + xlon[2] / 3600;
        float xlatdu = xlat[0] + xlat[1] / 60 + xlat[2] / 3600;
        float x = (Rlatdu < xlatdu) ? cal_distance(Rlatdu, Rlondu, xlatdu, Rlondu) : (-cal_distance(Rlatdu, Rlondu, xlatdu, Rlondu));
        float y = (Rlondu > xlondu) ? cal_distance(Rlatdu, Rlondu, Rlatdu, xlondu) : (-cal_distance(Rlatdu, Rlondu, Rlatdu, xlondu));
        float z = xalt - Ralt;
        // xyz2tpr
        float r = std::sqrt(x * x + y * y + z * z);
        float theta = 0;
        if ((x >= 0) && (y > 0)) theta = (x == 0) ? 90 : std::atan(y / x) * 180 / PIc;
        else if ((x >= 0) && (y < 0)) theta = (x == 0) ? -90 : std::atan(y / x) * 180 / PIc;
        else if ((x <= 0) && (y > 0)) theta = (x == 0) ? 90 : ((std::atan(y / x) * 180 / PIc) + 180);
        else if ((x <= 0) && (y < 0)) theta = (x == 0) ? -90 : ((std::atan(y / x) * 180 / PIc) - 180);
        if ((x >= 0) && (y == 0)) theta = 0;
        if ((x <= 0) && (y == 0)) theta = 180;
        float phi = (r == 0) ? 0 : std::asin(z / r) * 180 / PIc;
        // scan-range offsets
        float minradxtheta = 0, minradxphi = 0;
        float radxtheta = normalize_angle(theta), radxphi = phi;
        if (radxtheta < mintheta) minradxtheta = mintheta - radxtheta;
        else if (radxtheta > maxtheta) minradxtheta = radxtheta - maxtheta;
        if (radxphi < minphi) minradxphi = minphi - radxphi;
        else if (radxphi > maxphi) minradxphi = radxphi - maxphi;
        if (minradxtheta >= 90 || minradxphi >= 90) return 0;
        // PointPowerDensity (f = 24000 MHz -> choose_rDB returns 0 -> La = 1)
        float La = 1.0f;
        float etheta = tab[int(minradxtheta * PRECISION)];
        float ephi = tab[int(minradxphi * PRECISION)];
        float G_linear = std::pow(10, G / 10.0);
        float Gthetaphi = G_linear * std::abs(etheta) * std::abs(ephi);
        float Pt_watts = Pt * 1000.0;
        float Qt = (Pt_watts * Gthetaphi) / (4 * PIc * r * r * La);
        return Qt;
    }
};

struct CleanKernel {
    static constexpr double PIc = 3.14159265358979323846;
    static constexpr double EARTH_M = 6371193.0;
    Meta m; double Glin, PtW, dl;
    void init(const Meta& mm) { m = mm; Glin = std::pow(10.0, m.G / 10.0); PtW = m.Pt * 1000.0; dl = m.d / (300.0 / m.f); }
    double E(double offdeg) const {
        if (offdeg <= 0) return 1.0;
        double a = std::sin(offdeg * PIc / 180.0);
        double u = PIc * dl * a, den = m.N * std::sin(u);
        if (std::fabs(den) < 1e-300) return 1.0;
        return std::fabs(std::sin(m.N * u) / den);
    }
    double operator()(double lon, double lat, double alt) const {
        double rl1 = m.rlat * PIc / 180, rl2 = lat * PIc / 180;
        double x = EARTH_M * (rl2 - rl1);                                   // north
        double b = (m.rlon - lon) * PIc / 180;
        double y = EARTH_M * 2 * std::asin(std::cos(rl1) * std::sin(b / 2)); // west
        double z = alt - m.ralt;
        double r2 = x * x + y * y + z * z; if (r2 <= 0) return 0;
        double r = std::sqrt(r2);
        double theta = std::atan2(y, x) * 180 / PIc, phi = std::asin(z / r) * 180 / PIc;
        double dt = 0, dp = 0;
        if (theta < m.mintheta) dt = m.mintheta - theta; else if (theta > m.maxtheta) dt = theta - m.maxtheta;
        if (phi < m.minphi) dp = m.minphi - phi; else if (phi > m.maxphi) dp = phi - m.maxphi;
        if (dt >= 90 || dp >= 90) return 0;
        return PtW * Glin * E(dt) * E(dp) / (4 * PIc * r2);
    }
};
