set(SOURCES
    tests.cpp
)
add_executable(unit_tests ${SOURCES})

include(catch2)
target_link_libraries(unit_tests PUBLIC Catch2::Catch2WithMain)
target_include_directories(unit_tests PUBLIC ${PROJECT_SOURCE_DIR}/src)

set_property(TARGET unit_tests PROPERTY COMPILE_WARNING_AS_ERROR ON)

target_link_libraries(unit_tests PUBLIC bezier::bezier)
# target_include_directories(unit_tests PRIVATE ${IPRED_INCLUDE_DIR})
