#include <mpi.h>
#include "vendor/clipp.h"              //命令行参数解析库
#define CPPHTTPLIB_THREAD_POOL_COUNT 1 // single thread to avoid racial
#include "vendor/httplib.h"            //HTTP服务器库
#include "hiradar.hpp"                 //雷达类头文件
#include "radar_pool.hpp"              //引用雷达型号工厂，用于创建特定型号的雷达对象
#include "hiradar/TerrainGrid.hpp"
#include <gdal_priv.h>
#define PARAM_BUFFER_LENGTH 1024       // 定义一个宏，表示用于在MPI进程间传递参数的字符缓冲区大小
#ifndef RADAR
#define RADAR ANAPG66
#endif


// --- 2. handle函数：核心计算工作流 ---

// 这个函数代表了一次完整的计算任务
// 它接收一个雷达对象指针、剖分层级、探测范围、客户端传入的参数字符串，以及目标的RCS值
void handle(IRadar *radar, unsigned short level, float range, const char *params, float rcs)
{
    //std::cout << "[Rank ] Entered handle(). Parsing params..." << std::endl;
    Position pos;
    float phi = 0;
    float theta = 0;
    sscanf(params, "%f:%f:%f,%f:%f:%f,%f,%f,%f", pos.lon, pos.lon + 1, pos.lon + 2, pos.lat, pos.lat + 1, pos.lat + 2, &pos.alt, &phi, &theta);
    //建立以雷达为中心的网格（共100万个）
    auto point_list = CreateGrid(pos, range, level);
    auto pos_list = point_list.first;    // 获取指向网格点数组的指针
    auto pos_n = point_list.second;      // 获取网格点数量
    // 更新雷达对象的位置和姿态
    radar->Move(pos, phi, theta, 0);     
    // radar.PowerDensity(pos_list, pos_n);
    if (rcs > 0)
    {
        radar->CapableRegion(pos_list, pos_n, rcs);
    }
    else
    {
        //std::cout << "[Rank 1212] Calling CapablePowerDensity()..." << std::endl;
        radar->CapablePowerDensity(pos_list, pos_n);
    }
    delete[] pos_list;//释放内存
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    GDALAllRegister();
    int p_n, p_id;
    MPI_Comm_size(MPI_COMM_WORLD, &p_n);
    MPI_Comm_rank(MPI_COMM_WORLD, &p_id);

    caltools ct; //工具类对象
    Position center;

    std::string host = "0.0.0.0";
    int port = 10888;
    std::string db_host = "127.0.0.1";
    int db_port = 19000;
    std::string db_table = "radar";
    unsigned short level = 16;
    int range = 50000;
    int radar_id = 42;
    float rcs = -1; // if set rcs, calculate region instead

    auto cli =
        (clipp::option("-s") & clipp::value("rcs", rcs),
         clipp::option("-h", "--host") & clipp::value("host", host),
         clipp::option("-p", "--port") & clipp::value("port", port),
         clipp::option("-l", "--level") & clipp::value("level of the grid", level),
         clipp::option("-r", "--range") & clipp::value("maximum range to calculate", range),
         clipp::option("--db-host") & clipp::value("clickhouse host", db_host),
         clipp::option("--db-port") & clipp::value("clickhouse port", db_port),
         clipp::option("--db-table") & clipp::value("clickhouse table", db_table),
         clipp::value("radar id", radar_id));
    auto args = clipp::parse(argc, argv, cli);
    if (p_id == 0 && args.any_error())
    {
        for (auto &m : args.missing())
        {
            std::cerr << m.param()->label() << " not specified" << std::endl;
        }

        std::cerr << std::endl
                  << clipp::make_man_page(cli, argv[0]);
        exit(1);
    }
    // 使用了宏定义，调用雷达工厂函数，创建指定型号的雷达对象
    auto radar = RADAR(center, 0, 0);

    auto writer = new GeoSotChPowerDensityWriter(db_host, db_port, db_table, radar_id, level);
    radar->BindWriter(writer);
    TerrainGrid terrain_grid;
    bool success = terrain_grid.build("/mnt/d/ProgramData/document_keti/radar-demo1/final_cropped_area.tif");
    // bool success = false;
    if (!success)
    {
        std::cerr << "Failed to load terrain grid." << std::endl;
        exit(1);
    }
    radar->BindTerrain(&terrain_grid);
    char param_buf[PARAM_BUFFER_LENGTH] = {0};
    if (p_id == 0)
    {
        // 创建一个HTTP实例
        httplib::Server svr;

        // 定义当收到"/power"路径的PUT请求时的处理逻辑（回调函数）
        //  [...](){...}相当于一个lambda函数，
        svr.Put("/power", [&radar, &p_n, &range, &level, &param_buf, &rcs](const httplib::Request &req, httplib::Response &res)
                {
                    std::cout << "first callback" << std::endl;
                    //--- 以下是当请求到达时，主进程要执行的一系列动作 ---

                    // a. 获取任务参数
                    auto params = req.body.c_str(); // 从HTTP请求的正文中，获取客户端发送来的参数字符串
                    strcpy(param_buf, params);      // 将这些参数复制到一个用于MPI通信的共享缓冲区 param_buf 中

                    // b. 广播任务给所有从进程 (Workers)
                    // 循环遍历所有从进程（进程号从1到p_n-1）
                    for (size_t i = 1; i < p_n; i++)
                    {
                        MPI_Send(param_buf, PARAM_BUFFER_LENGTH, MPI_CHAR, i, 0, MPI_COMM_WORLD);
                    }
                    // c. 主进程自己也执行计算
                    // 调用handle函数，主进程亲自完成它自己那一份的计算任务
                    handle(radar, level, range, param_buf, rcs);

                    // d. 等待并回收所有从进程的完成信号
                    // 循环等待，直到接收到所有从进程发回的“完成”信号
                    for (size_t i = 1; i < p_n; i++)
                    {
                        int sig;
                        MPI_Status status;
                        // 一个进程调用此函数来等待并接收从其他进程发送过来的消息。
                        // 这是一个阻塞式 (Blocking) 操作，意味着如果还没有消息到达，程序会“卡”在这里，
                        // 直到收到消息才会继续往下执行。

                        // MPI_ANY_SOURCE 表示“我愿意接收来自任何一个从进程的消息”，而不在乎它们的顺序。
                        // PI_ANY_TAG 表示“我愿意接收任何主题的消息”。
                        // 1，MPI_INT代表要接收1个整形数据
                        // MPI_COMM_WORLD是一个通信域
                        MPI_Recv(&sig, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    }
                    std::cout << "Computation completed." << std::endl;
                });

        std::cout << "Radar Model started. Listening " << host << ":" << port << std::endl;
        // 启动服务器。这是一个阻塞式调用，程序会“卡”在这里，
        // 不会继续往下执行，直到有HTTP请求到来或程序被强制终止。
        // 这正是程序正常的工作状态：“待命”。
        svr.listen(host.c_str(), port);
        // if (!svr.listen(host.c_str(), port)) {
        //     std::cerr << "!!! FATAL ERROR: Failed to bind to port " << port << ". Is it already in use?" << std::endl;
        //     // 在这里添加MPI_Finalize()和退出逻辑
        //     MPI_Finalize();
        //     exit(1);
        // }
    }
    else
    {
        MPI_Status status;
        while (true)
        {
            MPI_Recv(param_buf, PARAM_BUFFER_LENGTH, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            handle(radar, level, range, param_buf, rcs);
            int sig = 0;
            MPI_Send(&sig, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
    }
    delete radar;
    delete writer;
    MPI_Finalize();
    return 0;
}