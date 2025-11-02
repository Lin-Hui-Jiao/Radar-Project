#include <random>
#include <iostream>
#include <fstream>
#include "geosot3d.hpp"
// #include <geosot3d/GeoSOT3D.hpp>
#include "vendor/clipp.h"

int main(int argc, char *argv[])
{
    //定义一个布尔变量，用于切换编码/解码模式。默认为false（编码模式）。
    bool reverse = false;

    int count = 1000;
    int level = 16;
    
    short code_type = 18;

    auto cli = (clipp::option("-r").set(reverse),
                clipp::option("-c") & clipp::value("count", count),
                clipp::option("-t") & clipp::value("code type", code_type),
                clipp::value("level", level));
    auto args = clipp::parse(argc, argv, cli);
    if (args.any_error())
    {
        for (auto &m : args.missing())
        {
            std::cerr << m.param()->label() << " not specified" << std::endl;
        }

        std::cerr << std::endl
                  << clipp::make_man_page(cli, argv[0]);
        exit(1);
    }
    // float lon = 0.01;
    // float lat = 0.01;
    // float alt = 2010;
    // std::cout << GeoSOT3D::Encode(lon, lat, alt, level) << std::endl;
    // std::cout << GeoSOT3D_Encode(lon, lat, alt, level) << std::endl;
    // return 0;

    srand(0);
    if (!reverse)
    {
        std::ofstream input;
        input.open("input.csv");
        std::ofstream output;
        output.open("output.csv");
        std::ofstream output2;
        output2.open("output2.csv");
        for (size_t i = 0; i < count; i++)
        {
            auto lon = rand() % 360 - 180;
            auto lat = rand() % 180 - 90;
            auto alt = (rand() % 50000000) / 1000.0;
            auto code = GeoSOT3D::Encode(lon, lat, alt, level);
            // auto code2 = GeoSOT3D_Encode(lon, lat, alt, level);
            input << i << "," << lon << "," << lat << "," << alt << "," << level << std::endl;
            output << i << "," << code << std::endl;
            // output2 << i << "," << code2 << std::endl;
        }
        input.close();
        output.close();
    }
    else
    {
        std::ofstream results;
        results.open("results.csv");
        std::ifstream input;
        input.open("output.csv");
        int id;
        float lon;
        float lat;
        float alt;
        char sep;
        for (size_t i = 0; i < count; i++)
        {
            uint64_t code;
            input >> id;
            input >> sep;
            input >> code;

            GeoSOT3D::Decode(lon, lat, alt, code, level);
            results << id << "," << lon << "," << lat << "," << alt << std::endl;
        }
        input.close();
        results.close();
    }

    return 0;
}
