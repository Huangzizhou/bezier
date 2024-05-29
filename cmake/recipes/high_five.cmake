if(TARGET HighFive)
    return()
endif()

message(STATUS "Third-party: creating target 'HighFive'")

include(CPM)
CPMAddPackage("gh:BlueBrain/HighFive#master")
