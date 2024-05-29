# Polyfem Solvers (https://github.com/polyfem/paraviewo)
# License: MIT

if(TARGET paraviewo::paraviewo)
    return()
endif()

message(STATUS "Third-party: creating target 'paraviewo::paraviewo'")

# ExternalProject_Add(paraviewo_external
#     GIT_REPOSITORY https://github.com/polyfem/paraviewo.git
#     GIT_TAG main
#     PREFIX ${CMAKE_CURRENT_BINARY_DIR}/paraviewo
#     CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/paraviewo
#     UPDATE_DISCONNECTED 1
# )

# set(PARAVIEWO_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/paraviewo/src/paraviewo_external/src)
# set(PARAVIEWO_LIBRARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/paraviewo/src/paraviewo_external-build)

# target_link_libraries(bezier PUBLIC ${PARAVIEWO_LIBRARY_DIR}/libparaviewo.a)
# target_include_directories(bezier PUBLIC ${PARAVIEWO_INCLUDE_DIR})

include(CPM)
CPMAddPackage("gh:polyfem/paraviewo#main")
