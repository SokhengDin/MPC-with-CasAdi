cmake_minimum_required(VERSION 3.9)
project(casadi)


find_package(PkgConfig)
pkg_search_module(CASADI casadi)

if(NOT CASADI_LIBRARIES)
  message(SEND_ERROR "Could not find casadi.pc - did you set PKG_CONFIG_PATH?")
endif()


message("CASADI_LIBRARY_DIRS: ${CASADI_LIBRARY_DIRS}")
message("CASADI_LIBRARIES: ${CASADI_LIBRARIES}")
message("CASADI_INCLUDEDIR: ${CASADI_INCLUDEDIR}")

# Add Ipopt link directories
if (WITH_IPOPT)
    link_directories(${IPOPT_LIBRARY_DIRS})
endif()

# Add executable
add_executable(test_linsol test_linsol.cpp)
target_link_libraries(test_linsol casadi)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")