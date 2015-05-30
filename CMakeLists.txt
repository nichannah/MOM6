cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(MOM6)

include_directories(${CMAKE_BINARY_DIR}/FMS config_src/dynamic src/framework)

file(GLOB SOURCES src/*/*.F90 src/*/*/*.F90 config_src/solo_driver/*.F90 src/framework/*.c)
add_executable(MOM6 ${SOURCES})
target_link_libraries(MOM6 fms)
target_link_libraries (MOM6 ${NETCDF_F90_LIBRARIES})
