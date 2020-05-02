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

add_library(st-sampler-lib
    random_spanning_trees.cpp
    approx_count_st.cpp
)

link_libraries(
    ${IGRAPH_LIBRARIES}
    st-sampler-lib
)

add_executable(test-igraph
    test-igraph.c
)

add_executable(st-sampler-full-graph-ratio main.cpp)
add_executable(st-sampler-sparse-graph-ratio main.cpp)
add_executable(st-sampler-ring-graph-ratio main.cpp)

add_executable(st-sampler-full-graph-variance main.cpp)
add_executable(st-sampler-sparse-graph-variance main.cpp)
add_executable(st-sampler-ring-graph-variance main.cpp)

add_executable(st-sampler-full-graph-constant main.cpp)
add_executable(st-sampler-sparse-graph-constant main.cpp)
add_executable(st-sampler-ring-graph-constant main.cpp)


# FULL GRAPH
target_compile_definitions(st-sampler-full-graph-ratio PRIVATE 
    FULL_GRAPH 
    CONVERGENCE_MODE=ApproxCountST::RATIO
)
target_compile_definitions(st-sampler-full-graph-variance PRIVATE 
    FULL_GRAPH 
    CONVERGENCE_MODE=ApproxCountST::VARIANCE
)
target_compile_definitions(st-sampler-full-graph-constant PRIVATE 
    FULL_GRAPH 
    CONVERGENCE_MODE=ApproxCountST::CONSTANT
)

# SPARSE GRAPH
target_compile_definitions(st-sampler-sparse-graph-ratio PRIVATE 
    SPARSE_GRAPH
    CONVERGENCE_MODE=ApproxCountST::RATIO
)
target_compile_definitions(st-sampler-sparse-graph-variance PRIVATE 
    SPARSE_GRAPH
    CONVERGENCE_MODE=ApproxCountST::VARIANCE
)
target_compile_definitions(st-sampler-sparse-graph-constant PRIVATE 
    SPARSE_GRAPH
    CONVERGENCE_MODE=ApproxCountST::CONSTANT
)

# RING GRAPH
target_compile_definitions(st-sampler-ring-graph-ratio PRIVATE 
    RING_GRAPH
    CONVERGENCE_MODE=ApproxCountST::RATIO
)
target_compile_definitions(st-sampler-ring-graph-variance PRIVATE 
    RING_GRAPH
    CONVERGENCE_MODE=ApproxCountST::VARIANCE
)
target_compile_definitions(st-sampler-ring-graph-constant PRIVATE 
    RING_GRAPH
    CONVERGENCE_MODE=ApproxCountST::CONSTANT
)

target_compile_definitions(st-sampler-lib PUBLIC
    PRESAMPLE_SIZE_REQUIRED=100
    PIVOT_BUFFER_SIZE=8
    RATIO_THRESHOLD_DEFAULT=0.002
    VARIANCE_THRESHOLD_DEFAULT=0.001
    CONSTANT_THRESHOLD_DEFAULT=4000
    INITIAL_REQUESTED_BATCH_SIZE=500
)

