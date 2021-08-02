CMAKE_MINIMUM_REQUIRED(VERSION 3.11)
PROJECT(XtraGraphComputing)
ADD_DEFINITIONS(-DUSE_LOG -DTIME_RESULTS)

FIND_PACKAGE(OpenMP REQUIRED)

if(NOT TARGET OpenMP::OpenMP_CXX)
    add_library(OpenMP_TARGET INTERFACE)
    add_library(OpenMP::OpenMP_CXX ALIAS OpenMP_TARGET random_walk/creeper.h)
    target_compile_options(OpenMP_TARGET INTERFACE ${OpenMP_CXX_FLAGS})
    find_package(Threads REQUIRED)
    target_link_libraries(OpenMP_TARGET INTERFACE Threads::Threads)
    target_link_libraries(OpenMP_TARGET INTERFACE ${OpenMP_CXX_FLAGS})
endif()

if(APPLE)
    SET(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -g -Xpreprocessor -fopenmp -lomp -march=native")
else()
    SET(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -g -fopenmp -march=native")
endif()


INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})

ADD_SUBDIRECTORY(util)
ADD_SUBDIRECTORY(toolset)
ADD_SUBDIRECTORY(random_walk)
