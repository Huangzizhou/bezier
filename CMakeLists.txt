# Detects whether this is a top-level project
get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if(HAS_PARENT)
    set(BEZIER_TOPLEVEL_PROJECT OFF)
else()
    set(BEZIER_TOPLEVEL_PROJECT ON)
endif()

cmake_minimum_required(VERSION 3.18)
project(bezier)
include(ExternalProject)


# Set the C++ standard
if(BEZIER_TOPLEVEL_PROJECT)
    set(CMAKE_CXX_STANDARD 20)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

# Options
option(IPRED_ARITHMETIC "Use the efficient Indirect Predicates library" ON)
option(UNIT_TESTS "Run unit tests" ON)
option(LAGVEC_GCC_O0 "Disable optimization for some complicated functions" ON)
option(HDF5_INTERFACE "Process HDF5 datasets" ON)
option(EIGEN_INTERFACE "Read data from Eigen classes" ON)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/recipes/")
include(bezier_cpm_cache)

if(LAGVEC_GCC_O0)
    add_compile_definitions(LAGVEC_GCC_O0)
endif()

# Set common compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall \
    -Wno-unknown-pragmas -fno-strict-aliasing -fopenmp -lpthread")

# Set release and debug specific compiler flags
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -pg -ggdb3")

add_library(bezier)
add_library(bezier::bezier ALIAS bezier)
set_target_properties(bezier PROPERTIES CXX_STANDARD 20)

# Use fast interval arithmetic
if(DEFINED IPRED_ARITHMETIC)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} \
        -O2 -frounding-math -msse2 -mavx2 -mfma")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} \
        -frounding-math -msse2 -mavx2 -mfma")
    set(USE_SIMD_INSTRUCTIONS ON)
else()
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
endif()

# Ensure the build type is set
if(NOT CMAKE_BUILD_TYPE)
  message(WARNING "CMAKE_BUILD_TYPE not specified. Defaulting to Release.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

##### EXECUTABLES ##############################################################
if (BEZIER_TOPLEVEL_PROJECT)
    add_executable(bin ${PROJECT_SOURCE_DIR}/src/main.cpp)
    target_link_libraries(bin PUBLIC bezier::bezier)
endif()

set(BEZIER_SOURCE_DIR "${PROJECT_SOURCE_DIR}/src")
set(BEZIER_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/src")

add_subdirectory(${BEZIER_SOURCE_DIR})

target_include_directories(bezier PUBLIC ${BEZIER_INCLUDE_DIR})

################################################################################
# Required libraries
################################################################################

##### Eigen ####################################################################
if (EIGEN_INTERFACE)
    add_compile_definitions(EIGEN_INTERFACE)
    include(eigen)
    target_link_libraries(bezier PUBLIC Eigen3::Eigen)
endif()

##### GMP ######################################################################
find_package(GMP REQUIRED)
target_link_libraries(bezier PUBLIC gmp::gmp)


##### INDIRECT_PREDICATES ######################################################
if (IPRED_ARITHMETIC)
    add_compile_definitions(IPRED_ARITHMETIC)

    include(indirect_predicates)
endif()


##### HDF5 #####################################################################
if (HDF5_INTERFACE)
    add_compile_definitions(HDF5_INTERFACE)
    find_package(HDF5 REQUIRED COMPONENTS C CXX)
    target_include_directories(bezier PUBLIC ${HDF5_INCLUDE_DIRS})
    target_link_libraries(bezier PUBLIC ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES})
endif()


################################################################################
# Tests
################################################################################

# Compile extras only if this is a top-level project
if(BEZIER_TOPLEVEL_PROJECT AND UNIT_TESTS)
    include(CTest)
    enable_testing()
    add_subdirectory(tests)
endif()