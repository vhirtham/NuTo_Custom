cmake_minimum_required(VERSION 3.9)
project(NuToApps)

FIND_PACKAGE(Boost REQUIRED unit_test_framework)
find_package(Eigen3)

find_path(EIGEN_INCLUDE_DIR NAMES Eigen/Core PATH_SUFFIXES eigen3)
include_directories(${EIGEN_INCLUDE_DIR})
include_directories("nuto/src")


include(ExternalProject)
ExternalProject_Add(NuTo
    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}"
               "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/nuto"
    BUILD_ALWAYS 1
    )


link_directories("${CMAKE_CURRENT_BINARY_DIR}/lib")

# ADD CUSTOM SCRIPTS
# ------------------
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/scripts/cmake ${CMAKE_MODULE_PATH})
include(IncludeScripts)

add_subdirectory(src)
add_subdirectory(tests)

