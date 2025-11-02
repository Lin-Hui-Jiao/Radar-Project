#include <string>
#include <iostream>
#include "vendor/httplib.h"
#include "vendor/clipp.h"

int main(int argc, char **argv)
{
    std::string host = "127.0.0.1";
    int port = 10888;
    std::string lon_str = "0:0:0";
    std::string lat_str = "0:0:0";
    float alt = 5000;
    float phi = 45;
    float theta = 30;

    auto cli =
        (clipp::option("-h", "--host") & clipp::value("host", host),
         clipp::option("-p", "--port") & clipp::value("port", port),
         clipp::value("longitude", lon_str),
         clipp::value("latitude", lat_str),
         clipp::value("altitude", alt),
         clipp::value("phi", phi),
         clipp::value("theta", theta));
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

    std::string conn_str;
    conn_str += "http://";
    conn_str += std::string(host);
    conn_str += ":";
    conn_str += std::to_string(port);
    std::cout << "Connecting radar server: " << conn_str << std::endl;
    // HTTP
    httplib::Client http(conn_str.c_str());

    std::stringstream ss;
    ss << lon_str << "," << lat_str << "," << alt << "," << phi << "," << theta;
    http.Put("/power", ss.str(), "text/plain");
    return 0;
}