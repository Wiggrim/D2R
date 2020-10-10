cmake_minimum_required(VERSION 2.8.3)
project(dop_rot)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g")

find_package(Eigen3)

include_directories( ${EIGEN3_INCLUDE_DIR} ${PROJECT_SOURCE_DIR} )

add_executable(dop_rot
    utility.cpp
    Filter.cpp
    UscentedFilter.cpp
    Dop2Rot.cpp
    main.cpp
)

target_link_libraries( dop_rot
    pthread
 )