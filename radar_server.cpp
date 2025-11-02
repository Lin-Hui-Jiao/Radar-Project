#include <mpi.h>
#include "vendor/clipp.h"
#define CPPHTTPLIB_THREAD_POOL_COUNT 1 // single thread to avoid racial
#include "vendor/httplib.h"
#include "hiradar/interface.hpp"
#include "hiradar/writer.hpp"
#include "hiradar/radar.hpp"
#include "hiradar/radar_pool.hpp"

#define PARAM_BUFFER_LENGTH 1024

void handle(IRadar *radar, unsigned short level, float range, const char *params)
{
    Position pos;
    float phi = 0;
    float theta = 0;
    sscanf(params, "%f:%f:%f,%f:%f:%f,%f,%f,%f", pos.lon, pos.lon + 1, pos.lon + 2, pos.lat, pos.lat + 1, pos.lat + 2, &pos.alt, &phi, &theta);
    //建立以雷达为中心的网格（共100万个）
    auto point_list = CreateGrid(pos, range, level);
    auto pos_list = point_list.first;
    auto pos_n = point_list.second;

    radar->Move(pos, phi, theta, 0);
    // radar.PowerDensity(pos_list, pos_n);
    // radar->CapablePowerDensity(pos_list, pos_n);
    radar->CapableRegion(pos_list, pos_n,7.5);
    free(pos_list);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int p_n, p_id;
    MPI_Comm_size(MPI_COMM_WORLD, &p_n);
    MPI_Comm_rank(MPI_COMM_WORLD, &p_id);

    caltools ct; //工具类对象
    Position center;

    std::string host = "0.0.0.0";
    int port = 10888;
    std::string db_host = "127.0.0.1";
    int db_port = 9000;
    std::string db_table = "radar";
    unsigned short level = 16;
    int range = 50000;
    int radar_id = 42;
    RadarParams params;

    auto cli =
        (clipp::option("-h", "--host") & clipp::value("host", host),
         clipp::option("-p", "--port") & clipp::value("port", port),
         clipp::option("-l", "--level") & clipp::value("level of the grid", level),
         clipp::option("-r", "--range") & clipp::value("maximum range to calculate", range),
         clipp::option("--db-host") & clipp::value("clickhouse host", db_host),
         clipp::option("--db-port") & clipp::value("clickhouse port", db_port),
         clipp::option("--db-table") & clipp::value("clickhouse table", db_table),
         clipp::value("阵元个数", params.N),
         clipp::value("阵元间距", params.d),
         clipp::value("工作频率", params.f),
         clipp::value("天线增益", params.G),
         clipp::value("峰值功率", params.Pt),
         clipp::value("最小方向角", params.mintheta),
         clipp::value("最大方向角", params.maxtheta),
         clipp::value("最小俯仰", params.minphi),
         clipp::value("最大俯仰", params.maxphi),
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

    auto radar = new Radar(params, center);

    auto writer = new GeoSotChPowerDensityWriter(db_host, db_port, db_table, radar_id, level);
    radar->BindWriter(writer);

    char param_buf[PARAM_BUFFER_LENGTH] = {0};
    if (p_id == 0)
    {
        // HTTP
        httplib::Server svr;
        svr.Put("/power", [&radar, &p_n, &range, &level, &param_buf](const httplib::Request &req, httplib::Response &res)
                {
                    auto params = req.body.c_str();
                    strcpy(param_buf, params);
                    for (size_t i = 1; i < p_n; i++)
                    {
                        MPI_Send(param_buf, PARAM_BUFFER_LENGTH, MPI_CHAR, i, 0, MPI_COMM_WORLD);
                    }
                    handle(radar, level, range, param_buf);
                    for (size_t i = 1; i < p_n; i++)
                    {
                        int sig;
                        MPI_Status status;
                        MPI_Recv(&sig, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    }
                    std::cout << "Computation completed." << std::endl;
                });

        std::cout << "Radar Model started. Listening " << host << ":" << port << std::endl;
        svr.listen(host.c_str(), port);
    }
    else
    {
        MPI_Status status;
        while (true)
        {
            MPI_Recv(param_buf, PARAM_BUFFER_LENGTH, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            handle(radar, level, range, param_buf);
            int sig = 0;
            MPI_Send(&sig, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
    }

    delete writer;
    MPI_Finalize();
    return 0;
}