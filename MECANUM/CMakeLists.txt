cmake_minimum_required(VERSION 3.11)
project(MECANUM CXX)
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)

find_package(PkgConfig)
find_package(casadi)
pkg_search_module(Eigen3 REQUIRED eigen3)
include_directories(${Eigen3_INCLUDE_DIRS})
pkg_search_module(casadi REQUIRED casadi)
include_directories(${casadi_INCLUDE_DIRS})
include_directories(include)

add_library(mpc_mecanum_lib src/auto_mecanum.cpp)
add_library(mpc_omni_lib src/auto_omni.cpp)

add_executable(mpc_controller main.cpp)
target_link_libraries(mpc_controller mpc_mecanum_lib casadi)

add_executable(mpc_omni main2.cpp)
target_link_libraries(mpc_omni mpc_omni_lib casadi)
