message(STATUS "Third-party: creating target 'zlib'")

include(CPM)
CPMAddPackage(
    NAME zlib_external
    URL https://www.zlib.net/current/zlib.tar.gz)
