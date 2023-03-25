cmake_minimum_required(VERSION 3.14)

project(MECANUM)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_BUILD_TYPE "Release")

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
find_package(Eigen3 REQUIRED)
find_package(casadi REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${Eigen3_INCLUDE_DIRS})
include_directories(include)


add_library(mpc_mecanum_lib src/mpc_mecanum.cpp)

add_executable(mpc_mecanum_controller main.cpp)
target_link_libraries(mpc_mecanum_controller mpc_mecanum_lib casadi)
