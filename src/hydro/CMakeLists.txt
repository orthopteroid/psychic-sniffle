cmake_minimum_required(VERSION 3.6)
project(hydro CXX)

file(GLOB LOCAL_SRC "*.cpp")

add_executable(hydro ${COMMON_SRC} ${LOCAL_SRC})
