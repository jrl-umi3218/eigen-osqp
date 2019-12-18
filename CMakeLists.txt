#
# Copyright 2018-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

cmake_minimum_required(VERSION 3.1.3)

set(PROJECT_NAME eigen-osqp)
set(PROJECT_DESCRIPTION "OSQP QP solver through Eigen3 library.")
set(PROJECT_URL "https://gite.lirmm.fr/multi-contact/eigen-osqp")
set(PROJECT_DEBUG_POSTFIX "_d")
set(PROJECT_VERSION 1.0.0)

set(CXX_DISABLE_WERROR True)
set(DOXYGEN_USE_MATHJAX "YES")

include(cmake/base.cmake)
include(cmake/boost.cmake)
include(cmake/msvc-specific.cmake)

option(COMPILE_EXAMPLES "Compile examples." ON)

project(${PROJECT_NAME} LANGUAGES CXX VERSION ${PROJECT_VERSION})

# For MSVC, set local environment variable to enable finding the built dll
# of the main library when launching ctest with RUN_TESTS
if(MSVC)
  set(CMAKE_MSVCIDE_RUN_PATH "\$(SolutionDir)/src/\$(Configuration)")
endif()

# add OSQP library
add_project_dependency(Eigen3 REQUIRED)
# FIXME We need OSQP >= 0.6.0 but this is not exported by their CMake files
add_project_dependency(osqp REQUIRED)

add_subdirectory(src)

if(${BUILD_TESTING})
  add_subdirectory(tests)
endif()

if(${COMPILE_EXAMPLES})
  add_subdirectory(examples)
endif()
