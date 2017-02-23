cmake_minimum_required(VERSION 3.6)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake-modules")

include_directories(inc)

add_subdirectory(src/schwefel)
add_subdirectory(src/hydro)
