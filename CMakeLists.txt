cmake_minimum_required(VERSION 3.6)
project(psychicsniffle)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake-modules")

include(GetCPUDetails)
get_cpu_details()

find_package(OpenMP)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build/${CMAKE_BUILD_TYPE}")

set(SOURCE_FILES taus88.h sniffle.h schwefel.h main.cpp)
add_executable(psychicsniffle ${SOURCE_FILES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if(CMAKE_BUILD_TYPE STREQUAL Debug)
	target_compile_definitions(psychicsniffle PRIVATE "-DDEBUG")

	target_compile_options(psychicsniffle PRIVATE "-g")
else()
	target_compile_options(psychicsniffle PRIVATE "-O3")

	if(OPENMP_FOUND)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif()

	if(CPU_HAS_SSE4)
		target_compile_options(psychicsniffle PRIVATE "-msse4")
	elseif(CPU_HAS_SSSE3)
		target_compile_options(psychicsniffle PRIVATE "-mssse3")
	elseif(CPU_HAS_SSE3)
		target_compile_options(psychicsniffle PRIVATE "-msse3")
	endif()
endif()
