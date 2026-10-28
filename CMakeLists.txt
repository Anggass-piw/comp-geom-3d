cmake_minimum_required(VERSION 3.16)
project(comp_geom_3d VERSION 0.1.0 LANGUAGES CXX)

option(GEOM3D_BUILD_TESTS "Build dependency-free tests" ON)

add_library(geom3d src/geometry.cpp)
target_include_directories(geom3d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(geom3d PUBLIC cxx_std_17)

if(MSVC)
    target_compile_options(geom3d PRIVATE /W4 /permissive-)
else()
    target_compile_options(geom3d PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(geom3d_demo app/main.cpp)
target_link_libraries(geom3d_demo PRIVATE geom3d)

if(GEOM3D_BUILD_TESTS)
    enable_testing()
    add_executable(geom3d_tests tests/test_geometry.cpp)
    target_link_libraries(geom3d_tests PRIVATE geom3d)
    add_test(NAME geom3d_tests COMMAND geom3d_tests)
endif()
