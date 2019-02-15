cmake_minimum_required(VERSION 2.8)

include(cmake/base.cmake)
include(cmake/boost.cmake)
include(cmake/eigen.cmake)
INCLUDE(cmake/msvc-specific.cmake)

set(PROJECT_NAME eigen-osqp)
set(PROJECT_DESCRIPTION "OSQP QP solver through Eigen3 library.")
set(PROJECT_URL "https://gite.lirmm.fr/multi-contact/eigen-osqp")

# Disable -Werror on Unix for now.
set(CXX_DISABLE_WERROR True)
set(DOXYGEN_USE_MATHJAX "YES")

option(USE_LTO "Enable link-time optimizations" ON)
CHECK_CXX_COMPILER_FLAG("-flto" CXX_HAS_LTO)
if(NOT CXX_HAS_LTO)
  set(USE_LTO OFF)
endif()
option(DISABLE_TESTS "Disable unit tests." OFF)
option(COMPILE_EXAMPLES "Compile examples." ON)

setup_project()

set(CXX_STANDARD 11)

if(UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
endif()

search_for_eigen()

# add OSQP library
find_package(osqp REQUIRED)

include_directories(src/)

add_subdirectory(src)

if(NOT ${DISABLE_TESTS})
  add_subdirectory(tests)
endif()

if(${COMPILE_EXAMPLES})
  add_subdirectory(examples)
endif()

# Add dependency towards the library in the pkg-config file.
pkg_config_append_libs(${PROJECT_NAME})
pkg_config_append_libs(osqp)

setup_project_finalize()
