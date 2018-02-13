cmake_minimum_required(VERSION 3.6)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake-modules")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_HOME_DIRECTORY}/${CMAKE_BUILD_TYPE}")

#######

include(GetCPUDetails)
get_cpu_details()

set(RELEASE_CXX_SSE "")
if(CPU_HAS_SSE4)
    set(RELEASE_CXX_SSE "-msse4")
elseif(CPU_HAS_SSSE3)
    set(RELEASE_CXX_SSE "-mssse3")
elseif(CPU_HAS_SSE3)
    set(RELEASE_CXX_SSE "-msse3")
endif()

#######

find_package(OpenMP)

########

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-multichar")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG ${RELEASE_CXX_SSE} ${OpenMP_CXX_FLAGS}")

#########

file(GLOB COMMON_SRC "src/*.cpp")

include_directories(src)

add_subdirectory(src/schwefel)
add_subdirectory(src/hydro)
add_subdirectory(src/quadratic)
