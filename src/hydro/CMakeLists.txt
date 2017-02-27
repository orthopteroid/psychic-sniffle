cmake_minimum_required(VERSION 3.6)
project(hydro CXX)

file(GLOB HYDRO_SRC "*.cpp" "../inc/hydro/*.h")

add_executable(hydro ${COMMON_SRC} ${HYDRO_SRC})
