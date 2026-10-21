# 雷达仿真系统 - 编译与执行指南

## 项目结构

```
radar-demo-0927-shiyan/
├── CMakeLists.txt          # 主构建配置文件
├── src/                    # 源代码目录
│   ├── caltools.cpp       # 计算工具实现
│   ├── radar.cpp          # 雷达核心功能
│   ├── occlusion_utils.cpp # 遮挡处理工具
│   ├── interface.cpp      # 接口实现
│   ├── server_new2.cpp    # 服务器主程序
│   ├── main.cpp           # 主程序入口
│   ├── test.cpp           # 测试程序
│   ├── Q-view.cpp         # 可视化程序
│   └── legacy/            # 旧版本代码
├── include/               # 头文件目录
├── geosot3d/             # GeoSOT3D 库
├── vendor/               # 第三方库
├── db/                   # 数据库脚本
└── build/                # 构建输出目录
    ├── bin/              # 可执行文件
    └── lib/              # 静态库文件
```

## 系统依赖

### 必需的软件包

```bash
# 基础编译环境
sudo apt-get update
sudo apt-get install build-essential cmake g++ gcc

# 线程和并行库
sudo apt-get install libpthread-stubs0-dev

# 地理信息处理库
sudo apt-get install libproj-dev libpng-dev

# ClickHouse 依赖
sudo apt-get install liblz4-dev libzstd-dev libcityhash-dev
```

### ClickHouse C++ 客户端（手动安装）

```bash
git clone https://github.com/ClickHouse/clickhouse-cpp.git
cd clickhouse-cpp
mkdir build && cd build
cmake ..
make -j$(nproc)
sudo make install
```

## 编译步骤

### 快速编译（推荐）

```bash
cd /home/rootjiao/radar-demo-0927-shiyan
rm -rf build
mkdir build && cd build
cmake ..
make -j4
```

### 详细步骤说明

**1. 清理旧构建**
```bash
rm -rf build
```

**2. 配置构建**
```bash
mkdir build && cd build
cmake ..
```

**3. 编译**
```bash
make -j4  # 使用4核并行编译
```

**4. 验证**
```bash
ls -lh bin/cityguardradar
```

## 运行程序

### 基本执行

```bash
cd build
./bin/cityguardradar -l 24 -r 300 -h 127.0.0.1 --db-port 19000 --db-table test.radar 42
```

### 参数说明

- `-l 24`: 层级（Level）
- `-r 300`: 半径（Radius，米）
- `-h 127.0.0.1`: 数据库主机
- `--db-port 19000`: 数据库端口
- `--db-table test.radar`: 数据库表名
- `--lon`: 经度（可选）
- `--lat`: 纬度（可选）
- `--height`: 雷达高度（可选，米）

### 执行示例

**指定经纬度：**
```bash
./bin/cityguardradar -l 24 -r 300 -h 127.0.0.1 --db-port 19000 \
  --lon 114.1669 --lat 22.2805 --db-table test.radar1 42
```

**指定高度：**
```bash
./bin/cityguardradar -l 25 -r 250 -h 127.0.0.1 --db-port 19000 \
  --lon 114.1670 --lat 22.2806 --height 30 --db-table test.radar7 42
```

**使用 GDB 调试：**
```bash
gdb --args ./bin/cityguardradar -l 24 -r 300 -h 127.0.0.1 \
  --db-port 19000 --db-table test.radar 42
```

## 常见问题

### 找不到 PROJ 库
```bash
sudo apt-get install libproj-dev
```

### 找不到 ClickHouse 库
从源码安装 clickhouse-cpp（见上文）

### 编译内存不足
```bash
make -j2  # 减少并行数
```

### 运行时找不到动态库
```bash
sudo ldconfig
# 或
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
```

## 构建产物

- `build/bin/cityguardradar` - 主程序
- `build/lib/libhiradar.a` - 核心库
- `build/lib/libgeosot3d.a` - GeoSOT3D库

## 开发提示

- 修改源代码后直接运行 `make`
- 修改 CMakeLists.txt 后需重新运行 `cmake ..`
- 调试模式已启用（CMAKE_BUILD_TYPE=Debug）
