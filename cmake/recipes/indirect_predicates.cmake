message(STATUS "Third-party: creating target 'indirect_predicates'")

include(CPM)
CPMAddPackage(
    NAME indirect_predicates_external
    GIT_REPOSITORY https://github.com/fsichetti/Indirect_Predicates.git
    GIT_TAG master)

message(STATUS "indirect_predicates_external_SOURCE_DIR: ${indirect_predicates_external_SOURCE_DIR}")
target_include_directories(bezier PUBLIC ${indirect_predicates_external_SOURCE_DIR}/include/)
