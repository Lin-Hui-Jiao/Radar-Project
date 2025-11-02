# ==============================================================================
# 1. Basic Project Setup
# ==============================================================================

cmake_minimum_required(VERSION 3.10)
project(RadarSimulation CXX)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set output directories for build artifacts
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# ==============================================================================
# 2. Find External Dependencies
# ==============================================================================

# --- Standard Packages ---
find_package(Threads REQUIRED) # For pthread
find_package(OpenMP REQUIRED) # <-- 【修改】从 MPI 改为 OpenMP
find_package(PNG REQUIRED) # <-- 【新增】寻找 libpng 库
# --- 【修改】手动寻找 PROJ 库 ---
# 如果 find_package 失败，我们手动来找
find_path(PROJ_INCLUDE_DIR NAMES proj.h HINTS /usr/include /usr/local/include)
find_library(PROJ_LIBRARY NAMES proj HINTS /usr/lib/x86_64-linux-gnu /usr/local/lib)

if(NOT PROJ_INCLUDE_DIR OR NOT PROJ_LIBRARY)
    message(FATAL_ERROR "Could not find PROJ library headers or binaries. Please run 'sudo apt-get install libproj-dev'")
else()
    # 创建一个导入的目标，方便后续使用
    add_library(PROJ::proj UNKNOWN IMPORTED)
    set_target_properties(PROJ::proj PROPERTIES
        IMPORTED_LOCATION "${PROJ_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${PROJ_INCLUDE_DIR}")
    message(STATUS "Found PROJ: ${PROJ_LIBRARY}")
endif()

# --- Manually Find clickhouse-cpp and its dependencies (like in the Makefile) ---
# This section replaces the failing `find_package(clickhouse-cpp)`
find_library(CLICKHOUSE_CPP_LIBRARY NAMES clickhouse-cpp-lib HINTS /usr/local/lib /usr/lib)
find_library(LZ4_LIBRARY NAMES lz4 HINTS /usr/local/lib /usr/lib)
find_library(ZSTD_LIBRARY NAMES zstd HINTS /usr/local/lib /usr/lib)
find_library(CITYHASH_LIBRARY NAMES cityhash HINTS /usr/local/lib /usr/lib)

# Check if the manual search was successful
if(NOT CLICKHOUSE_CPP_LIBRARY)
    message(FATAL_ERROR "Could not find libclickhouse-cpp-lib.a or .so. Please check installation.")
endif()

# # 【核心修改】将所有源文件集中到一个列表中，不再创建中间库
# set(ALL_SOURCES
#     server_new2.cpp
#     radar.cpp
#     caltools.cpp
#     occlusion_utils.cpp
# )

# add_executable(radar ${ALL_SOURCES})
# # target_compile_definitions(city_radar_sim PRIVATE RADAR=CityGuardRadar)

# # 为可执行文件配置头文件搜索路径
# target_include_directories(radar PUBLIC
#     ${CMAKE_PROJECT_SOURCE_DIR}/include
#     ${CMAKE_PROJECT_SOURCE_DIR}/vendor  # 存放 httplib.h, clipp.h 等
# )

# # 为可执行文件链接所有需要的库
# target_link_libraries(radar PUBLIC
#     # 并行库
#     OpenMP::OpenMP_CXX   # <-- 【修改】链接 OpenMP
#     Threads::Threads
#     geosot3d
#     # 地理信息库
#     PROJ::proj           # <-- 【修改】确保链接 PROJ

#     # 数据库库
#     ${CLICKHOUSE_CPP_LIBRARY}
#     ${LZ4_LIBRARY}
#     ${ZSTD_LIBRARY}
#     ${CITYHASH_LIBRARY}
# )
# install(TARGETS radar DESTINATION bin)

# ==============================================================================
# 3. Add Subdirectory for geosot3d Library
# ==============================================================================

# Build the geosot3d library from its own subdirectory
# This assumes 'geosot3d/CMakeLists.txt' exists and creates the 'geosot3d' target
 add_subdirectory(geosot3d)

# ==============================================================================
# 4. Build the Core 'hiradar' Library
# ==============================================================================

# Define the source files for the hiradar library
set(HIRADAR_SRCS
    src/caltools.cpp
    src/radar.cpp
    src/occlusion_utils.cpp
    src/interface.cpp
)

# Create the hiradar static library target
add_library(hiradar STATIC ${HIRADAR_SRCS})

# Configure include paths for the hiradar library
target_include_directories(hiradar PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/vendor
)

# Link hiradar against the geosot3d library
# 声明 hiradar 库需要链接 geosot3d 库，并且它自身也需要MPI和GDAL
target_link_libraries(hiradar PUBLIC
    geosot3d
    PROJ::proj
    OpenMP::OpenMP_CXX
    PNG::PNG
)

# ==============================================================================
# 5. Build the Final Server Executables
# ==============================================================================

# List of radar models to generate executables for
set(RADAR_MODELS CityGuardRadar)
# ANAPS115 ANAPY12 ANAPY2 ANSPS49 SPY1D ANAPS145 ANMPQ65，这个地方是有问题的
# Loop through each model to create a specific server executable
foreach(MODEL ${RADAR_MODELS})
    string(TOLOWER "${MODEL}" TARGET_NAME)

    add_executable(${TARGET_NAME} src/server_new2.cpp)

    # Add the -DRADAR=... preprocessor definition for this target
    target_compile_definitions(${TARGET_NAME} PRIVATE RADAR=${MODEL})

    # Link the executable against all necessary libraries
    target_link_libraries(${TARGET_NAME} PRIVATE
    # 只需要链接 hiradar，CMake会自动处理其 PUBLIC 依赖
    hiradar
    # --- 只链接那些 server.cpp 直接使用，但 hiradar 没有传递的库 ---  
    # 假设 server.cpp 直接使用了 clickhouse 和 httplib (需要线程)
    ${CLICKHOUSE_CPP_LIBRARY}
    ${LZ4_LIBRARY}
    ${ZSTD_LIBRARY}
    ${CITYHASH_LIBRARY}
    Threads::Threads
    OpenMP::OpenMP_CXX
    PNG::PNG
)
endforeach()

# ==============================================================================
# 6. Build DEM Comparison Experiment Module
# ==============================================================================

# Find GDAL library for DEM generation
find_package(GDAL REQUIRED)

# Define source files for comparison experiment
# TODO: comparison_experiment is temporarily disabled due to missing VisibilityComparator class
# set(COMPARISON_SRCS
#     src/dem_generator.cpp
#     src/dem_loader.cpp
#     src/dem_visibility.cpp
#     src/comparison_experiment.cpp
# )

# Create comparison experiment executable
# add_executable(comparison_experiment ${COMPARISON_SRCS})

# Configure include paths
# target_include_directories(comparison_experiment PRIVATE
#     ${CMAKE_CURRENT_SOURCE_DIR}/include
#     ${CMAKE_CURRENT_SOURCE_DIR}/vendor
#     ${GDAL_INCLUDE_DIRS}
# )

# Link against all necessary libraries
# target_link_libraries(comparison_experiment PRIVATE
#     hiradar              # Core radar library (includes geosot3d, PROJ, OpenMP)
#     ${GDAL_LIBRARIES}    # GDAL for DEM generation
#     PNG::PNG             # PNG output (optional, for visualization)
# )

# Install comparison experiment executable
# install(TARGETS comparison_experiment DESTINATION bin)

# message(STATUS "Configured comparison_experiment executable")
# message(STATUS "  GDAL include: ${GDAL_INCLUDE_DIRS}")
# message(STATUS "  GDAL libraries: ${GDAL_LIBRARIES}")

# ==============================================================================
# 7. Build Fixed-Height Comparison Experiment Module
# ==============================================================================

# Define source files for fixed-height experiment
set(FIXED_HEIGHT_SRCS
    src/fixed_height_comparator.cpp
    src/dem_generator.cpp
    src/dem_loader.cpp
    src/dem_visibility.cpp
    src/fixed_height_experiment.cpp
)

# Create fixed-height experiment executable
add_executable(fixed_height_experiment ${FIXED_HEIGHT_SRCS})

# Configure include paths
target_include_directories(fixed_height_experiment PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/vendor
    ${GDAL_INCLUDE_DIRS}
)

# Link against all necessary libraries
target_link_libraries(fixed_height_experiment PRIVATE
    hiradar              # Core radar library (includes geosot3d, PROJ, OpenMP)
    ${GDAL_LIBRARIES}    # GDAL for DEM generation and loading
    PNG::PNG             # PNG output for visualization
)

# Install fixed-height experiment executable
install(TARGETS fixed_height_experiment DESTINATION bin)

message(STATUS "Configured fixed_height_experiment executable")
message(STATUS "  GDAL include: ${GDAL_INCLUDE_DIRS}")
message(STATUS "  GDAL libraries: ${GDAL_LIBRARIES}")