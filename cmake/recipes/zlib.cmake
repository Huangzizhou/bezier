# ExternalProject_Add(zlib_external
#     URL https://www.zlib.net/current/zlib.tar.gz
#     PREFIX ${CMAKE_CURRENT_BINARY_DIR}/zlib
#     CONFIGURE_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/zlib/src/zlib_external/configure --prefix=${CMAKE_CURRENT_BINARY_DIR}/zlib/install --static
#     BUILD_COMMAND make
#     BUILD_IN_SOURCE 1
#     INSTALL_COMMAND make install
#     UPDATE_DISCONNECTED 1
# )

# set(ZLIB_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/zlib/install/include)
# set(ZLIB_LIBRARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/zlib/install/lib)

message(STATUS "Third-party: creating target 'zlib'")

include(CPM)
CPMAddPackage(
    NAME zlib_external
    URL https://www.zlib.net/current/zlib.tar.gz)
