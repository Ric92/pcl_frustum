cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(frustum)

set (CMAKE_CXX_STANDARD 11)

find_package(PCL 1.2 REQUIRED)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})
link_directories(${EIGEN3_LIBRARY_DIRS})
add_definitions(${EIGEN3_DEFINITIONS})

add_executable (frustum main.cpp)
target_link_libraries (frustum ${PCL_LIBRARIES})
target_link_libraries (frustum ${EIGEN3_LIBRARIES})
