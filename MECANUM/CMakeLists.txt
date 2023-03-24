cmake_minimum_required(VERSION 3.14)
project(MECANUM)

include_directories(include)

find_package(casadi REQUIRED)
find_package(eigen3 REQUIRED)

add_executable(mpc_mecanum src/mpc_mecanum)
target_link_libraries(mpc_mecanum casadi)