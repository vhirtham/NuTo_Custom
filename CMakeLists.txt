cmake_minimum_required(VERSION 3.7)
project(NuToApps)
set(CMAKE_VERBOSE_MAKEFILE OFF)

option(ENABLE_TESTS true)

enable_testing()

set(NUTO_APPS_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++14 -Wall -Wno-long-long -pedantic")

FIND_PACKAGE(Boost REQUIRED unit_test_framework)

find_path(EIGEN_INCLUDE_DIR NAMES Eigen/Core PATH_SUFFIXES eigen3)
include_directories(${EIGEN_INCLUDE_DIR})
include_directories("nuto")


include(ExternalProject)
ExternalProject_Add(NuTo_ext
    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}"
               "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
               #"-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/nuto"
    BUILD_ALWAYS 1
    )


link_directories("${CMAKE_CURRENT_BINARY_DIR}/lib")

# ADD CUSTOM SCRIPTS
# ------------------
set(CMAKE_MODULE_PATH ${NUTO_APPS_ROOT_DIR}/scripts/cmake ${CMAKE_MODULE_PATH})
set(CMAKE_MODULE_PATH ${NUTO_APPS_ROOT_DIR}/nuto/scripts/cmake ${CMAKE_MODULE_PATH})
include(IncludeScripts)

add_subdirectory(applications)
add_subdirectory(src)
add_subdirectory(tests)
