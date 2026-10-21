# field_export — standalone build (does not touch the repository's top-level CMakeLists.txt).
#   cmake -S field_export -B field_export/build -DCMAKE_BUILD_TYPE=Release && cmake --build field_export/build -j
cmake_minimum_required(VERSION 3.16)
project(field_export CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
# Same floating-point behaviour as the prototypes that were checked bit-exact against the CSVs: no FMA contraction,
# no fast-math.
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
add_compile_options(-ffp-contract=off)

find_package(OpenMP REQUIRED)
find_library(PROJ_LIB proj REQUIRED)

set(HIRADAR_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/..)
add_library(fe_occlusion STATIC ${HIRADAR_ROOT}/src/occlusion_utils.cpp)
target_include_directories(fe_occlusion PUBLIC ${HIRADAR_ROOT}/include ${HIRADAR_ROOT}/include/hiradar)

foreach(t field_export field_check)
  add_executable(${t} ${t}.cpp)
  target_link_libraries(${t} PRIVATE fe_occlusion OpenMP::OpenMP_CXX ${PROJ_LIB})
  target_compile_definitions(${t} PRIVATE FIELD_EXPORT_SRC_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
endforeach()
