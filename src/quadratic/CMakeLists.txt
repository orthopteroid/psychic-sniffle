cmake_minimum_required(VERSION 3.6)
project(quadratic CXX)

file(GLOB LOCAL_SRC "*.cpp")

add_executable(quadratic ${COMMON_SRC} ${LOCAL_SRC})
