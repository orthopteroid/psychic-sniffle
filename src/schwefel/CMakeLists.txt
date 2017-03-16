cmake_minimum_required(VERSION 3.6)
project(schwefel CXX)

file(GLOB LOCAL_SRC "*.cpp")

add_executable(schwefel ${COMMON_SRC} ${LOCAL_SRC})
