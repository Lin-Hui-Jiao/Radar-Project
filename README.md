# 安装

- 安装基础编译环境 gcc,pthread 等
- 安装 [clickhouse-cpp](https://github.com/ClickHouse/clickhouse-cpp)，[absl](https://github.com/abseil/abseil-cpp)
- 安装 [mpich](https://www.mpich.org/)
- 安装 geosot3d(随代码提供)

```bash
mkdir -p build
make -j 8  # 编译库，简单的 server/client 封装
sudo make install  # 安装库与头文件到 /usr/local 下
```

八型雷达：

先按上一步安装库，然后进入radar目录：
```bash
mkdir -p build
make -j 8
```

# 运行

运行八型雷达（radar目录下）
```bash
bash db/create_db.sh http://127.0.0.1:18123 test.radar && mpiexec -n 12 ./cityguardradar -l 25 -r 800 -h 127.0.0.1 --db-port 19000 --db-table test.radar 42
```

 ./cityguardradar -l 24 -r 300 -h 127.0.0.1 --db-port 19000 --db-table test.radar 42
 gdb --args ./cityguardradar -l 24 -r 300 -h 127.0.0.1 --db-port 19000 --db-table test.radar 42


./cityguardradar -l 24 -r 300 -h 127.0.0.1 --db-port 19000 --lon 114.1669 --lat 22.2805 --db-table test.radar1 42


./cityguardradar -l 24 -r 300 -h 127.0.0.1 --db-port 19000 --lon 114.1667 --lat 22.2810 --db-table test.radar3 42


 ./cityguardradar -l 25 -r 250 -h 127.0.0.1 --db-port 19000 --lon 114.1670 --lat 22.2806 --db-table test.radar7 --height 30 42