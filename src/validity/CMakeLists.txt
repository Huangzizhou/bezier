set(HEADERS
    element_validity.hpp
    lagrangeVector.hpp
    transMatrices.hpp
)

set(SOURCES
    element_validity.cpp
    lagrangeVector.cpp
    transMatrices.cpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(bezier PUBLIC ${SOURCES})
