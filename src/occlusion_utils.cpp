#include "hiradar/occlusion_utils.h"

bool lineSegIntersectTri(Vec3d line[2], Vec3d tri[3], Vec3d *point)
{
    // Calculate the two edges of the triangle and the direction vector of the line segment.
    Vec3d e0 = tri[1] - tri[0];
    Vec3d e1 = tri[2] - tri[0];
    Vec3d dir = line[1] - line[0];
    Vec3d dir_norm = dir.normalize();
    Vec3d h = dir_norm.cross(e1);
    const double a = e0.dot(h);

    if (a > -EPSILON && a < EPSILON)
    {
        return false;
    }

    Vec3d s = line[0] - tri[0];
    const double f = 1.0f / a;
    const double u = f * s.dot(h);

    if (u < 0.0f || u > 1.0f)
    {
        return false;
    }

    Vec3d q = s.cross(e0);
    const double v = f * dir_norm.dot(q);

    if (v < 0.0f || u + v > 1.0f)
    {
        return false;
    }

    const double t = f * e1.dot(q);

    if (t > EPSILON && t < sqrtf(dir.dot(dir))) // segment intersection
    {
        if (point)
        {
            *point = line[0] + dir_norm * t;
        }

        return true;
    }

    return false;
}


// 这是配合 R-tree 索引库的回调函数，用于确定给定 (x, y) 位置的高度。
// This function is a callback used to determine the height at a given (x, y) location.
double HeightCallback(ValueType id, const double x, const double y)
{
    // Create a triangle from the given id, which contains vertex positions.
    Vec3d tri[3] =
        {
            {get<0>(id), get<1>(id), get<2>(id)},
            {get<3>(id), get<4>(id), get<5>(id)},
            {get<6>(id), get<7>(id), get<8>(id)},
        };
    // std::cout << "id = ("
    //           << std::get<0>(id) << ", " << std::get<1>(id) << ", " << std::get<2>(id) << ", "
    //           << std::get<3>(id) << ", " << std::get<4>(id) << ", " << std::get<5>(id) << ", "
    //           << std::get<6>(id) << ", " << std::get<7>(id) << ", " << std::get<8>(id)
    //           << ")" << std::endl;

    // Define a line segment from (x, y, MAXH) to (x, y, MINH).
    Vec3d line[2] =
        {
            {x, y, MAXH},
            {x, y, MINH},
    };
    Vec3d *point = new Vec3d();
    // Check if the line segment intersects with the triangle.
    if (lineSegIntersectTri(line, tri, point))
    {
        return point->z;
    }
    // If there is no intersection, return a value indicating to continue searching.
    return MINH - 1;
}

// This function is a callback used to check for intersection between a line segment and a triangle.
// 这是配合 R-tree 索引库的回调函数，用于检查给定线段是否与三角形相交。只要有一个相交，直接返回false，就不用进一步搜索了。
bool IntersectCallback(ValueType id, const double a[3], const double b[3])
{
    Vec3d tri[3] =
        {
            {get<0>(id), get<1>(id), get<2>(id)},
            {get<3>(id), get<4>(id), get<5>(id)},
            {get<6>(id), get<7>(id), get<8>(id)},
        };
    Vec3d line[2] =
        {
            {a[0], a[1], a[2]},
            {b[0], b[1], b[2]},
        };
    // If an intersection is found, return false to indicate the intersection.
    if (lineSegIntersectTri(line, tri, NULL))
    {
        return false;
    }
    // If no intersection is found, return true to indicate continue searching.
    return true;
}

void Usage()
{
    printf("Usage:   [--index:     input file path  ]\n");
}
float DMSToDecimal(const float dms[3])
{
    return dms[0] + dms[1] / 60.0f + dms[2] / 3600.0f;
}







