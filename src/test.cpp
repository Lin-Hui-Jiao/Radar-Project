#include <iostream>
#include <cmath>     // For std::abs
#include <iomanip>   // For std::fixed and std::setprecision
#include <proj.h>    // PROJ library header

// 定义一个函数来比较两个浮点数是否在误差范围内相等
bool areAlmostEqual(double a, double b, double epsilon = 1e-9) {
    return std::abs(a - b) < epsilon;
}

int main() {
    // --- 1. 初始化：定义所有常量和输入坐标 ---
    std::cout << std::fixed << std::setprecision(10); // 设置输出精度

    // 原始 WGS84 经纬度坐标 (EPSG:4326)
    const double initial_lon_deg = 114.1670;
    const double initial_lat_deg = 22.2806;

    // 坐标系常量，与您项目中的完全一致
    const double min_x = 835000.0;
    const double min_y = 815500.0;

    std::cout << "--- 初始坐标 (WGS84) ---" << std::endl;
    std::cout << "Longitude: " << initial_lon_deg << std::endl;
    std::cout << "Latitude:  " << initial_lat_deg << std::endl;
    std::cout << "\n" << std::endl;

    // --- 2. 执行正向转换 (WGS84 -> 局部坐标) ---
    std::cout << "--- 步骤 1: 正向转换 (WGS84 -> 局部坐标) ---" << std::endl;

    double local_x = 0.0;
    double local_y = 0.0;

    PJ_CONTEXT *ctx_fwd = proj_context_create();
    PJ *proj_fwd = proj_create_crs_to_crs(ctx_fwd, "EPSG:4326", "EPSG:2326", NULL);
    if (proj_fwd == nullptr) {
        std::cerr << "正向投影创建失败!" << std::endl;
        return 1;
    }

    PJ_COORD wgs84_in = proj_coord(initial_lat_deg, initial_lon_deg, 0, 0);
    PJ_COORD hk1980_out = proj_trans(proj_fwd, PJ_FWD, wgs84_in);

    // 注意：PROJ输出的xy顺序与我们通常理解的相反
    // hk1980_out.xy.x 是 Northing (北向值), 对应纬度
    // hk1980_out.xy.y 是 Easting (东向值), 对应经度
    local_x = hk1980_out.xy.y - min_x;
    local_y = hk1980_out.xy.x - min_y;
    
    std::cout << "转换后的香港格网坐标 (Easting, Northing): " 
              << hk1980_out.xy.y << ", " << hk1980_out.xy.x << std::endl;
    std::cout << "计算出的局部坐标 (local_x, local_y): " 
              << local_x << ", " << local_y << std::endl;
    std::cout << "\n" << std::endl;

    proj_destroy(proj_fwd);
    proj_context_destroy(ctx_fwd);

    // --- 3. 执行逆向转换 (局部坐标 -> WGS84) ---
    std::cout << "--- 步骤 2: 逆向转换 (局部坐标 -> WGS84) ---" << std::endl;
    
    double final_lon_deg = 0.0;
    double final_lat_deg = 0.0;

    PJ_CONTEXT *ctx_inv = proj_context_create();
    PJ *proj_inv = proj_create_crs_to_crs(ctx_inv, "EPSG:2326", "EPSG:4326", NULL);
    if (proj_inv == nullptr) {
        std::cerr << "逆向投影创建失败!" << std::endl;
        return 1;
    }

    // 严格的逆过程：先恢复偏移，得到完整的香港格网坐标
    double restored_easting = local_x + min_x;
    double restored_northing = local_y + min_y;
    std::cout << "从局部坐标恢复的香港格网坐标 (Easting, Northing): "
              << restored_easting << ", " << restored_northing << std::endl;

    PJ_COORD hk1980_in = proj_coord(restored_northing, restored_easting, 0, 0);
    PJ_COORD wgs84_final_out = proj_trans(proj_inv, PJ_FWD, hk1980_in);

    // 对于WGS84(EPSG:4326)，输出顺序是 (Latitude, Longitude)
    final_lat_deg = wgs84_final_out.xy.x;
    final_lon_deg = wgs84_final_out.xy.y;

    std::cout << "最终转换回的 WGS84 坐标:" << std::endl;
    std::cout << "Longitude: " << final_lon_deg << std::endl;
    std::cout << "Latitude:  " << final_lat_deg << std::endl;
    std::cout << "\n" << std::endl;

    proj_destroy(proj_inv);
    proj_context_destroy(ctx_inv);

    // --- 4. 验证结果 ---
    std::cout << "--- 步骤 3: 结果验证 ---" << std::endl;
    if (areAlmostEqual(initial_lon_deg, final_lon_deg) && areAlmostEqual(initial_lat_deg, final_lat_deg)) {
        std::cout << "✅ 验证成功！正向和逆向过程互为逆运算。" << std::endl;
    } else {
        std::cout << "❌ 验证失败！两个过程不是严格的逆运算。" << std::endl;
        std::cout << "经度差: " << std::abs(initial_lon_deg - final_lon_deg) << std::endl;
        std::cout << "纬度差: " << std::abs(initial_lat_deg - final_lat_deg) << std::endl;
    }

    return 0;
}