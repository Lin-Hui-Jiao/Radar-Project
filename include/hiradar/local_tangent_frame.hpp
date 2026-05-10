#pragma once

#include "hiradar/occlusion_utils.h"

namespace hiradar {

class LocalTangentFrame {
public:
    LocalTangentFrame();
    LocalTangentFrame(double origin_lon_deg, double origin_lat_deg, double origin_h_m);

    void Reset(double origin_lon_deg, double origin_lat_deg, double origin_h_m);

    bool IsInitialized() const;

    Vec3d ToENU(double lon_deg, double lat_deg, double h_m) const;
    Vec3d FromENU(double east_m, double north_m, double up_m) const;

    double origin_lon_deg() const { return origin_lon_deg_; }
    double origin_lat_deg() const { return origin_lat_deg_; }
    double origin_h_m() const { return origin_h_m_; }

private:
    static Vec3d GeodeticToECEF(double lon_deg, double lat_deg, double h_m);
    static Vec3d ECEFToGeodetic(const Vec3d& ecef);

    double origin_lon_deg_;
    double origin_lat_deg_;
    double origin_h_m_;
    double sin_lat_;
    double cos_lat_;
    double sin_lon_;
    double cos_lon_;
    Vec3d origin_ecef_;
    bool initialized_;
};

}  // namespace hiradar
