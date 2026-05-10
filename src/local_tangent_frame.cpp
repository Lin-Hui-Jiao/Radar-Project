#include "hiradar/local_tangent_frame.hpp"

#include <cmath>
#include <stdexcept>

namespace hiradar {

namespace {

constexpr double kWgs84A = 6378137.0;
constexpr double kWgs84F = 1.0 / 298.257223563;
constexpr double kWgs84B = kWgs84A * (1.0 - kWgs84F);
constexpr double kWgs84E2 = kWgs84F * (2.0 - kWgs84F);
constexpr double kWgs84Ep2 = (kWgs84A * kWgs84A - kWgs84B * kWgs84B) / (kWgs84B * kWgs84B);
constexpr double kPi = 3.14159265358979323846;

double DegToRad(double degrees) {
    return degrees * kPi / 180.0;
}

double RadToDeg(double radians) {
    return radians * 180.0 / kPi;
}

}  // namespace

LocalTangentFrame::LocalTangentFrame()
    : origin_lon_deg_(0.0),
      origin_lat_deg_(0.0),
      origin_h_m_(0.0),
      sin_lat_(0.0),
      cos_lat_(1.0),
      sin_lon_(0.0),
      cos_lon_(1.0),
      origin_ecef_(0.0, 0.0, 0.0),
      initialized_(false) {}

LocalTangentFrame::LocalTangentFrame(double origin_lon_deg, double origin_lat_deg, double origin_h_m)
    : LocalTangentFrame() {
    Reset(origin_lon_deg, origin_lat_deg, origin_h_m);
}

void LocalTangentFrame::Reset(double origin_lon_deg, double origin_lat_deg, double origin_h_m) {
    origin_lon_deg_ = origin_lon_deg;
    origin_lat_deg_ = origin_lat_deg;
    origin_h_m_ = origin_h_m;

    const double lat_rad = DegToRad(origin_lat_deg_);
    const double lon_rad = DegToRad(origin_lon_deg_);
    sin_lat_ = std::sin(lat_rad);
    cos_lat_ = std::cos(lat_rad);
    sin_lon_ = std::sin(lon_rad);
    cos_lon_ = std::cos(lon_rad);
    origin_ecef_ = GeodeticToECEF(origin_lon_deg_, origin_lat_deg_, origin_h_m_);
    initialized_ = true;
}

bool LocalTangentFrame::IsInitialized() const {
    return initialized_;
}

Vec3d LocalTangentFrame::ToENU(double lon_deg, double lat_deg, double h_m) const {
    if (!initialized_) {
        throw std::runtime_error("LocalTangentFrame is not initialized");
    }

    const Vec3d ecef = GeodeticToECEF(lon_deg, lat_deg, h_m);
    const double dx = ecef.x - origin_ecef_.x;
    const double dy = ecef.y - origin_ecef_.y;
    const double dz = ecef.z - origin_ecef_.z;

    const double east = -sin_lon_ * dx + cos_lon_ * dy;
    const double north = -sin_lat_ * cos_lon_ * dx - sin_lat_ * sin_lon_ * dy + cos_lat_ * dz;
    const double up = cos_lat_ * cos_lon_ * dx + cos_lat_ * sin_lon_ * dy + sin_lat_ * dz;
    return Vec3d(east, north, up);
}

Vec3d LocalTangentFrame::FromENU(double east_m, double north_m, double up_m) const {
    if (!initialized_) {
        throw std::runtime_error("LocalTangentFrame is not initialized");
    }

    const double dx = -sin_lon_ * east_m - sin_lat_ * cos_lon_ * north_m + cos_lat_ * cos_lon_ * up_m;
    const double dy = cos_lon_ * east_m - sin_lat_ * sin_lon_ * north_m + cos_lat_ * sin_lon_ * up_m;
    const double dz = cos_lat_ * north_m + sin_lat_ * up_m;
    const Vec3d ecef(origin_ecef_.x + dx, origin_ecef_.y + dy, origin_ecef_.z + dz);
    return ECEFToGeodetic(ecef);
}

Vec3d LocalTangentFrame::GeodeticToECEF(double lon_deg, double lat_deg, double h_m) {
    const double lon_rad = DegToRad(lon_deg);
    const double lat_rad = DegToRad(lat_deg);
    const double sin_lat = std::sin(lat_rad);
    const double cos_lat = std::cos(lat_rad);
    const double sin_lon = std::sin(lon_rad);
    const double cos_lon = std::cos(lon_rad);
    const double prime_vertical = kWgs84A / std::sqrt(1.0 - kWgs84E2 * sin_lat * sin_lat);

    const double x = (prime_vertical + h_m) * cos_lat * cos_lon;
    const double y = (prime_vertical + h_m) * cos_lat * sin_lon;
    const double z = (prime_vertical * (1.0 - kWgs84E2) + h_m) * sin_lat;
    return Vec3d(x, y, z);
}

Vec3d LocalTangentFrame::ECEFToGeodetic(const Vec3d& ecef) {
    const double p = std::sqrt(ecef.x * ecef.x + ecef.y * ecef.y);
    const double lon = std::atan2(ecef.y, ecef.x);
    const double theta = std::atan2(ecef.z * kWgs84A, p * kWgs84B);
    const double sin_theta = std::sin(theta);
    const double cos_theta = std::cos(theta);

    const double lat = std::atan2(
        ecef.z + kWgs84Ep2 * kWgs84B * sin_theta * sin_theta * sin_theta,
        p - kWgs84E2 * kWgs84A * cos_theta * cos_theta * cos_theta
    );

    const double sin_lat = std::sin(lat);
    const double prime_vertical = kWgs84A / std::sqrt(1.0 - kWgs84E2 * sin_lat * sin_lat);
    const double h = p / std::cos(lat) - prime_vertical;

    return Vec3d(RadToDeg(lon), RadToDeg(lat), h);
}

}  // namespace hiradar
