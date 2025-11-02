#include "httplib.h"            //HTTP服务器库
#include "clipp.h"              //命令行参数解析库
//#include "hiradar.hpp"                 //雷达类头文件
#include "hiradar/writer.hpp"
#include "hiradar/radar_pool.hpp"              //引用雷达型号工厂，用于创建特定型号的雷达对象
// #include <gdal_priv.h>
#include "hiradar/grid.hpp"
#ifndef RADAR
#define RADAR CityGuardRadar
#endif


int main(int argc, char **argv)
{
    std::string host = "0.0.0.0";
    int port = 10888;
    std::string db_host = "127.0.0.1";
    int db_port = 19000;
    std::string db_table = "radar";
    int range = 250;
    int radar_id = 42;
    int level = 25;
    float rcs = -1; // if set rcs, calculate region instead
    float lon = 114.1670;
    float lat = 22.2806;
    float alt = 80;
    float fixed_height = 10.0; // 默认固定高度为10米
    auto cli =
        (clipp::option("-s") & clipp::value("rcs", rcs),
         clipp::option("-h", "--host") & clipp::value("host", host),
         clipp::option("-p", "--port") & clipp::value("port", port),
         clipp::option("-l", "--level") & clipp::value("level of the grid", level),
         clipp::option("-r", "--range") & clipp::value("maximum range to calculate", range),
         clipp::option("--height") & clipp::value("fixed calculation height in meters", fixed_height),
         clipp::option("--db-host") & clipp::value("clickhouse host", db_host),
         clipp::option("--db-port") & clipp::value("clickhouse port", db_port),
         clipp::option("--db-table") & clipp::value("clickhouse table", db_table),
         clipp::option("--lon") & clipp::value("lon", lon),
         clipp::option("--lat") & clipp::value("lat", lat),
         clipp::option("--alt") & clipp::value("alt", alt),
         clipp::value("radar id", radar_id));
    auto args = clipp::parse(argc, argv, cli);
    Position radar_pos;
    DecimalToDMS(radar_pos.lon, lon);
    DecimalToDMS(radar_pos.lat, lat);
    radar_pos.alt = alt;

    float base_phi = 0.0f;   // 载具基础俯仰角
    float base_theta = 0.0f; // 载具基础方位角
    // 使用了宏定义，调用雷达工厂函数，创建指定型号的雷达对象
    auto radar = RADAR(radar_pos, base_phi, base_theta);

    auto writer = new GeoSotChPowerDensityWriter(db_host, db_port, db_table, radar_id, level);
    radar->BindWriter(writer);
    RTree3d* rtree = new RTree3d();
    if (!rtree->Load("/mnt/d/ProgramData/document_keti/radar-demo-0916/test_area.3idx")) {
        std::cerr << "Fatal Error: Failed to load R-tree index file." << std::endl;
        return 1;
    }
    cout << "三角面片已经加载完毕！"<<endl;
    cout << "range: " << range << ", level: " << level << ", fixed_height: " << fixed_height << "m" << endl;
    radar->BindRtree(rtree);
    cout << "雷达经度是:" << lon << " " << radar_pos.lon[0] << " " << radar_pos.lon[1] << " " << radar_pos.lon[2] << endl;
    cout << "雷达的纬度是" << lat << endl;
    radar->Move(radar_pos, base_phi, base_theta, 0);
    if (rcs > 0) {
        //radar->CapableRegion(level, range, rcs);
        cout << "目前暂不支持该功能，待后续开发！" << endl;
        return 1;
    } else {
        // radar->CapablePowerDensity(range,level);
        // radar->CapableGroundPowerDensity(range,1.0);
        for(int i = 1;i < 30;i += 2){
            radar->CapableFixedHeightPowerDensity(range, 1.0, i);
        }
        
    }
    std::cout << "--> All tasks completed." << std::endl;

    // --- 4. 清理资源 ---
    std::cout << "--> Cleaning up resources..." << std::endl;
    delete radar;
    delete writer;
    delete rtree;
    return 0;
}