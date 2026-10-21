# vehicle_field512 — standalone build of the checking tool (does not touch the repository's top-level CMakeLists.txt).
# The fields themselves are computed by ../field_export (built separately, see README.md).
#   cmake -S vehicle_field512 -B vehicle_field512/build -DCMAKE_BUILD_TYPE=Release && cmake --build vehicle_field512/build -j
cmake_minimum_required(VERSION 3.16)
project(vehicle_field512 CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
add_compile_options(-ffp-contract=off)

find_package(OpenMP REQUIRED)
find_library(PROJ_LIB proj REQUIRED)

set(HIRADAR_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/..)
add_library(vf_occlusion STATIC ${HIRADAR_ROOT}/src/occlusion_utils.cpp)
target_include_directories(vf_occlusion PUBLIC ${HIRADAR_ROOT}/include ${HIRADAR_ROOT}/include/hiradar)

add_executable(site_check site_check.cpp)
target_link_libraries(site_check PRIVATE vf_occlusion OpenMP::OpenMP_CXX ${PROJ_LIB})

# independent visibility check (RTree::Search + watertight segment / triangle test), see vis_check.cpp
add_executable(vis_check vis_check.cpp)
target_include_directories(vis_check PRIVATE ${HIRADAR_ROOT}/field_export)
target_link_libraries(vis_check PRIVATE vf_occlusion OpenMP::OpenMP_CXX ${PROJ_LIB})
