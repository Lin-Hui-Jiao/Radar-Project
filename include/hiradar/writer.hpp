#pragma once

#include <clickhouse/client.h>
#include <memory>
#include <geosot3d.hpp>
#include "hiradar/interface.hpp"
// #include "hiradar/grid.hpp"
#include "hiradar/occlusion_utils.h"

namespace ch = clickhouse;


// 定义一个专门用于将GeoSOT编码和功率密度数据写入ClickHouse的类
// 它公开继承自 IPowerDensityWriter 接口
class GeoSotChPowerDensityWriter : public IPowerDensityWriter
{
public://table是数据库表格的名字
    GeoSotChPowerDensityWriter(const std::string &host, int port, const std::string &table, uint64_t radar_id, unsigned short level)
        : _client(new ch::Client(ch::ClientOptions().SetHost(host).SetPort(port))), _table(table), _radar_id(radar_id), _level(level) {}
        // C++成员初始化列表：在函数体执行前，初始化所有成员变量
        // 关键：在这里根据传入的host和port，创建了一个ClickHouse客户端连接对象 _client

     // 虚函数Write的实现，这是该类的核心工作函数    
    virtual void Write(const Position *pos_list, float *values, size_t n, float timestamp)
    {

        // 1. 准备数据块
        // Block 是 clickhouse-cpp 库中用于批量插入数据的容器
        ch::Block block;

         // 2. 创建并填充每一列的数据
        // 为 "id" 列创建一个64位无符号整型列对象
        auto id = std::make_shared<ch::ColumnUInt64>();
        for (size_t i = 0; i < n; i++)
        {
            id->Append(_radar_id); // 将同一个 radar_id 添加 n 次
        }
        // 为 "timestamp" 列创建一个32位浮点数列对象
        auto ts = std::make_shared<ch::ColumnFloat32>();
        for (size_t i = 0; i < n; i++)
        {
            
            ts->Append(timestamp);
        }
        // 创建经度、纬度、高程和编码的列对象

        auto lon = std::make_shared<ch::ColumnFloat32>();
        auto lat = std::make_shared<ch::ColumnFloat32>();
        auto alt = std::make_shared<ch::ColumnFloat32>();
        auto code = std::make_shared<ch::ColumnUInt64>();
        // 遍历输入的n个位置点数据
        for (size_t i = 0; i < n; i++)
        {
            auto pos = pos_list[i];
            auto lon_val = DMSToDecimal(pos.lon);
            auto lat_val = DMSToDecimal(pos.lat);
            code->Append(GeoSOT3D::Encode(pos.lon, pos.lat, pos.alt, _level)); //这个地方可能编码还是有点问题的！
            lon->Append(lon_val);
            lat->Append(lat_val);
            alt->Append(pos.alt);
        }

        auto value = std::make_shared<ch::ColumnFloat32>();
        for (size_t i = 0; i < n; i++)
        {
            value->Append(values[i]);
        }
        // 3. 将所有填充好的列组装进数据块(Block)中
        // 参数一是数据库中的列名，参数二是对应的列数据对象

        block.AppendColumn("id", id);
        block.AppendColumn("code", code);
        block.AppendColumn("lon", lon);
        block.AppendColumn("lat", lat);
        block.AppendColumn("alt", alt);
        block.AppendColumn("density", value);
        block.AppendColumn("timestamp", ts);
        // 4. 执行数据库插入操作
        // 调用客户端的Insert方法，将整个数据块一次性地写入到指定的表中
        _client->Insert(_table, block);
    }

    void SetLevel(unsigned short level) { _level = level; }

private:
    std::unique_ptr<ch::Client> _client;
    std::string _table;
    uint64_t _radar_id;
    unsigned short _level;
};
