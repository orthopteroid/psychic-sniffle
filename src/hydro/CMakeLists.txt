cmake_minimum_required(VERSION 3.6)
project(hydro CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake-modules")

include(GetCPUDetails)
get_cpu_details()

find_package(OpenMP)
find_library(LIBRT_LIBRARIES rt)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories(../inc)
file(GLOB COMMON_SRC "../src/*.cpp" "../inc/*.h")

file(GLOB HYDRO_SRC "*.cpp" "../inc/hydro/*.h")

add_executable(hydro ${COMMON_SRC} ${HYDRO_SRC})

target_link_libraries(hydro "-lrt")

if(CMAKE_BUILD_TYPE STREQUAL Debug)
	target_compile_definitions(hydro PRIVATE "-DDEBUG")

	target_compile_options(hydro PRIVATE "-g")
else()
	target_compile_options(hydro PRIVATE "-O3")

	if(OPENMP_FOUND)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif()

	if(CPU_HAS_SSE4)
		target_compile_options(hydro PRIVATE "-msse4")
	elseif(CPU_HAS_SSSE3)
		target_compile_options(hydro PRIVATE "-mssse3")
	elseif(CPU_HAS_SSE3)
		target_compile_options(hydro PRIVATE "-msse3")
	endif()
endif()

