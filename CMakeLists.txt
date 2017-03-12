cmake_minimum_required(VERSION 3.3)
project(tom)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(SOURCE_FILES src/tom.cpp include/tom)

include_directories(include/external)
include_directories(include/tom)
include_directories(/opt/local/include/eigen3)

add_custom_target(deploy COMMAND make -C ${tom_SOURCE_DIR} deploy)
add_custom_target(deploy_mkl_openmp COMMAND make -C ${tom_SOURCE_DIR} deploy_mkl_openmp)

add_executable(tom ${SOURCE_FILES})
