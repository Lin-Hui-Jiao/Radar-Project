#ifndef OCCLUSION_UTILS_H
#define OCCLUSION_UTILS_H

#include "query3dRtree.h"
#include "kernel.h"
#include <proj.h> //地理空间转换库
using namespace std;
typedef tuple<double, double, double, double, double, double, double, double, double> ValueType;
typedef RTree<ValueType, double, 3> RTree3d;

struct Rect3d
{
    Rect3d() {}

    Rect3d(double minX, double minY, double minZ, double maxX, double maxY, double maxZ)
    {
        min[0] = minX;
        min[1] = minY;
        min[2] = minZ;
        max[0] = maxX;
        max[1] = maxY;
        max[2] = maxZ;
    }
    double min[3];
    double max[3];
};
// various operations
// 非常方便的进行向量运算
class Vec3d
{
public:
    double x, y, z;
    // Calculate the dot product of this vector with another vector b.
    double dot(const Vec3d &b)
    {
        return Vec3d::x * b.x + Vec3d::y * b.y + Vec3d::z * b.z;
    }
    // Calculate the cross product of this vector with another vector b.
    Vec3d cross(const Vec3d &b)
    {
        return Vec3d(
            Vec3d::y * b.z - Vec3d::z * b.y,
            Vec3d::z * b.x - Vec3d::x * b.z,
            Vec3d::x * b.y - Vec3d::y * b.x);
    }
    // Normalize this vector to make its length equal to 1.
    Vec3d normalize()
    {
        const double s = 1.0f / sqrtf(Vec3d::x * Vec3d::x + Vec3d::y * Vec3d::y + Vec3d::z * Vec3d::z);
        return Vec3d(Vec3d::x * s, Vec3d::y * s, Vec3d::z * s);
    }
    Vec3d operator+(const Vec3d &b)
    {
        return Vec3d(
            Vec3d::x + b.x,
            Vec3d::y + b.y,
            Vec3d::z + b.z);
    }
    Vec3d operator+=(const Vec3d &b)
    {
        *this = Vec3d::operator+(b);
        return *this;
    }
    Vec3d operator-(const Vec3d &b)
    {
        return Vec3d(
            Vec3d::x - b.x,
            Vec3d::y - b.y,
            Vec3d::z - b.z);
    }
    Vec3d operator-=(const Vec3d &b)
    {
        *this = Vec3d::operator-(b);
        return *this;
    }
    Vec3d operator*(const Vec3d &b)
    {
        return Vec3d(
            Vec3d::x * b.x,
            Vec3d::y * b.y,
            Vec3d::z * b.z);
    }
    Vec3d operator*=(const Vec3d &b)
    {
        *this = Vec3d::operator*(b);
        return *this;
    }
    Vec3d operator*(double b)
    {
        return Vec3d(
            Vec3d::x * b,
            Vec3d::y * b,
            Vec3d::z * b);
    }
    Vec3d operator*=(double b)
    {
        *this = Vec3d::operator*(b);
        return *this;
    }
    Vec3d operator/(const Vec3d &b)
    {
        return Vec3d(
            Vec3d::x / b.x,
            Vec3d::y / b.y,
            Vec3d::z / b.z);
    }
    Vec3d operator/=(const Vec3d &b)
    {
        *this = Vec3d::operator/(b);
        return *this;
    }
    Vec3d operator/(double b)
    {
        return Vec3d(
            Vec3d::x * b,
            Vec3d::y * b,
            Vec3d::z * b);
    }
    Vec3d operator/=(double b)
    {
        *this = Vec3d::operator/(b);
        return *this;
    }
    Vec3d(double x, double y, double z)
    {
        Vec3d::x = x;
        Vec3d::y = y;
        Vec3d::z = z;
    }
    Vec3d(double x)
    {
        Vec3d::x = x;
        Vec3d::y = x;
        Vec3d::z = x;
    }
    Vec3d()
    {
        //
    }
    ~Vec3d()
    {
        //
    }
};
// 一个简单的类，用于管理PROJ资源的生命周期 (RAII模式)
struct PROJ_Manager {
    PJ_CONTEXT *ctx = nullptr;
    PJ *proj_from = nullptr;

    PROJ_Manager() {
        ctx = proj_context_create();
        proj_from = proj_create_crs_to_crs(ctx, "EPSG:4326", "EPSG:2326", NULL);
    }

    ~PROJ_Manager() {
        if (proj_from) proj_destroy(proj_from);
        if (ctx) proj_context_destroy(ctx);
    }
};

#define EPSILON 0.000001f
const double PI = 3.14159265358979323846;
bool lineSegIntersectTri(Vec3d line[2], Vec3d tri[3], Vec3d *point);
double HeightCallback(ValueType id, const double x, const double y);
bool IntersectCallback(ValueType id, const double a[3], const double b[3]);
void Usage();
float DMSToDecimal(const float dms[3]);









#endif // OCCLUSION_UTILS_H
