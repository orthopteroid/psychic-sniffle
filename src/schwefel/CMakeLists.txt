cmake_minimum_required(VERSION 3.6)
project(schwefel CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake-modules")

include(GetCPUDetails)
get_cpu_details()

find_package(OpenMP)
find_library(LIBRT_LIBRARIES rt)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories(../inc)
file(GLOB COMMON_SRC "../src/*.cpp" "../inc/*.h")

file(GLOB SCHWEFEL_SRC "*.cpp" "../inc/schwefel/*.cpp")

add_executable(schwefel ${COMMON_SRC} ${SCHWEFEL_SRC})

if(CMAKE_BUILD_TYPE STREQUAL Debug)
	target_compile_definitions(schwefel PRIVATE "-DDEBUG")

	target_compile_options(schwefel PRIVATE "-g")
else()
	target_compile_options(schwefel PRIVATE "-O3")

	if(OPENMP_FOUND)
		#		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif()

	if(CPU_HAS_SSE4)
		target_compile_options(schwefel PRIVATE "-msse4")
	elseif(CPU_HAS_SSSE3)
		target_compile_options(schwefel PRIVATE "-mssse3")
	elseif(CPU_HAS_SSE3)
		target_compile_options(schwefel PRIVATE "-msse3")
	endif()
endif()
