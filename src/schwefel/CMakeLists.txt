cmake_minimum_required(VERSION 3.6)
project(schwefel CXX)

file(GLOB SCHWEFEL_SRC "*.cpp" "../inc/schwefel/*.cpp")

add_executable(schwefel ${COMMON_SRC} ${SCHWEFEL_SRC})
