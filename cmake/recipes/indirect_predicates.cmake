# ExternalProject_Add(indirect_predicates_external
#     GIT_REPOSITORY https://github.com/MarcoAttene/Indirect_Predicates.git
#     GIT_TAG master
#     PREFIX ${CMAKE_CURRENT_BINARY_DIR}/indirect_predicates
#     CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/indirect_predicates
#     UPDATE_DISCONNECTED 1
#     INSTALL_COMMAND ""
#     TEST_COMMAND ""
# )

# set(IPRED_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/indirect_predicates/src/indirect_predicates_external/include)
# target_include_directories(bezier PRIVATE ${IPRED_INCLUDE_DIR})

message(STATUS "Third-party: creating target 'indirect_predicates'")

include(CPM)
CPMAddPackage(
    NAME indirect_predicates_external
    GIT_REPOSITORY https://github.com/MarcoAttene/Indirect_Predicates.git
    GIT_TAG master)

target_include_directories(bezier PUBLIC ${indirect_predicates_external_SOURCE_DIR}/include/)
