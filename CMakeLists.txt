cmake_minimum_required(VERSION 3.10)
project(st-stampler-igraph)

# enable C++17 standard (current directory scope)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# configure build types and build c-flags
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
# set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3") # without -DNDEBUG

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/")
find_package(IGRAPH REQUIRED)

find_package(Eigen3 REQUIRED)

include_directories(
    include
    ${IGRAPH_INCLUDES}
    ${EIGEN3_INCLUDE_DIRS}
)

link_libraries(
    ${IGRAPH_LIBRARIES}
)

add_executable(test-igraph
    test-igraph.c
)

add_executable(st-sampler-full-graph
    random_spanning_trees.cpp
    approx_count_st.cpp
    main.cpp
)

add_executable(st-sampler-sparse-graph
    random_spanning_trees.cpp
    approx_count_st.cpp
    main.cpp
)

add_executable(st-sampler-ring-graph
    random_spanning_trees.cpp
    approx_count_st.cpp
    main.cpp
)

target_compile_definitions(st-sampler-full-graph PUBLIC FULL_GRAPH)
target_compile_definitions(st-sampler-sparse-graph PUBLIC SPARSE_GRAPH)
target_compile_definitions(st-sampler-ring-graph PUBLIC RING_GRAPH)


