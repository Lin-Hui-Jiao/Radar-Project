# 将geosot3d.cpp 源文件编译成一个名为 "geosot3d" 的静态库
add_library(geosot3d STATIC geosot3d.cpp)

# 因为geosot3d.cpp 源文件中使用了其他库，所以需要将其他库也编译成静态库
target_include_directories(geosot3d PUBLIC
    # ${CMAKE_CURRENT_SOURCE_DIR} 是一个CMake变量，代表当前CMakeLists.txt所在的目录
    ${CMAKE_CURRENT_SOURCE_DIR}         # 用于找到 geosot3d.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/vendor  # 用于找到 libmorton
)
