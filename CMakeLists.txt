cmake_minimum_required(VERSION 3.1.3)

set(PROJECT_NAME eigen-osqp)
set(PROJECT_DESCRIPTION "OSQP QP solver through Eigen3 library.")
set(PROJECT_URL "https://gite.lirmm.fr/multi-contact/eigen-osqp")
set(PROJECT_DEBUG_POSTFIX "_d")
set(PROJECT_VERSION 0.4.0)

include(cmake/base.cmake)
include(cmake/boost.cmake)
include(cmake/msvc-specific.cmake)

# Disable -Werror on Unix for now.
# set(CXX_DISABLE_WERROR True)
set(DOXYGEN_USE_MATHJAX "YES")

option(USE_LTO "Enable link-time optimizations" ON)
CHECK_CXX_COMPILER_FLAG("-flto" CXX_HAS_LTO)
if(NOT CXX_HAS_LTO)
  set(USE_LTO OFF)
endif()
option(COMPILE_EXAMPLES "Compile examples." ON)

project(${PROJECT_NAME} LANGUAGES CXX VERSION ${PROJECT_VERSION})

# For MSVC, set local environment variable to enable finding the built dll
# of the main library when launching ctest with RUN_TESTS
IF(MSVC)
  SET(CMAKE_MSVCIDE_RUN_PATH "\$(SolutionDir)/src/\$(Configuration)")
ENDIF(MSVC)

if(UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
endif()

# add OSQP library
add_project_dependency(Eigen3 REQUIRED)
add_project_dependency(osqp REQUIRED ">= 0.6.0")

add_subdirectory(src)

if(${BUILD_TESTING})
  add_subdirectory(tests)
endif()

if(${COMPILE_EXAMPLES})
  add_subdirectory(examples)
endif()

# Add dependency towards the library in the pkg-config file.
pkg_config_append_libs(${PROJECT_NAME})
pkg_config_append_libs(osqp)

setup_project_finalize()
setup_project_package_finalize()
