set(HEADERS
    element_validity.hpp
)

set(SOURCES
    element_validity.cpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(bezier PUBLIC ${SOURCES})
