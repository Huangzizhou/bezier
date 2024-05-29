# Detects whether this is a top-level project
get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if(HAS_PARENT)
    set(JACOBIAN_TOPLEVEL_PROJECT OFF)
else()
    set(JACOBIAN_TOPLEVEL_PROJECT ON)
endif()

cmake_minimum_required(VERSION 3.18)
project(jacobian)
include(ExternalProject)

# Options
option(IPRED_ARITHMETIC "Use the efficient Indirect Predicates library" ON)
option(JACOBIAN_WITH_TESTS "Run unit tests" OFF)

set(CMAKE_CXX_FLAGS "-Wall -Wno-unknown-pragmas -fno-strict-aliasing -fopenmp -lpthread")
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg -ggdb3")

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/recipes/")

include(jacobian_cpm_cache)

# Use fast interval arithmetic
if(IPRED_ARITHMETIC)
    set(CMAKE_CXX_FLAGS_RELEASE "-frounding-math -O2 -msse2 -mavx2")
    set(USE_SIMD_INSTRUCTIONS ON)
else()
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()

# Set the C++ standard
if(JACOBIAN_TOPLEVEL_PROJECT)
    set(CMAKE_CXX_STANDARD 20)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()


add_library(jacobian)
add_library(jacobian::jacobian ALIAS jacobian)

set(JACOBIAN_SOURCE_DIR "${PROJECT_SOURCE_DIR}/src")
set(JACOBIAN_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/src")

# Add source and header files to jacobian
add_subdirectory(${JACOBIAN_SOURCE_DIR})

# Include folders
target_include_directories(jacobian PUBLIC ${JACOBIAN_INCLUDE_DIR})

################################################################################
# Required libraries
################################################################################

##### Eigen ####################################################################
# include(eigen)
# target_link_libraries(jacobian PUBLIC Eigen3::Eigen)

##### GMP ######################################################################
find_package(GMP REQUIRED)
target_link_libraries(jacobian PUBLIC gmp::gmp)


##### INDIRECT_PREDICATES ######################################################
if (IPRED_ARITHMETIC)
    add_compile_definitions(IPRED_ARITHMETIC)

    include(indirect_predicates)
endif()

if (JACOBIAN_TOPLEVEL_PROJECT)
    add_executable(jacobian_bin ${PROJECT_SOURCE_DIR}/src/main.cpp)
    target_link_libraries(jacobian_bin PUBLIC jacobian::jacobian)

    # add_executable(jacobian_test_bin ${PROJECT_SOURCE_DIR}/src/test.cpp)
    # target_link_libraries(jacobian_test_bin PUBLIC jacobian::jacobian)
endif()


################################################################################
# Tests
################################################################################

# Compile extras only if this is a top-level project
if(JACOBIAN_TOPLEVEL_PROJECT AND JACOBIAN_WITH_TESTS)
    include(CTest)
    enable_testing()
    add_subdirectory(tests)
endif()