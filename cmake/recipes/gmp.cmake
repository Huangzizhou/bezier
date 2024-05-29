ExternalProject_Add(gmp_external
    URL https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gmp
    CONFIGURE_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/gmp/src/gmp_external/configure --prefix=${CMAKE_CURRENT_BINARY_DIR}/gmp/install --enable-cxx > /dev/null
    BUILD_COMMAND make
    BUILD_IN_SOURCE 1
)

set(GMP_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/gmp/install/include)
set(GMP_LIBRARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/gmp/install/lib)
target_link_libraries(jacobian PRIVATE ${GMP_LIBRARY_DIR}/libgmp.a ${GMP_LIBRARY_DIR}/libgmpxx.a)
target_include_directories(jacobian PRIVATE ${GMP_INCLUDE_DIR})
