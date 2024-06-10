set(HEADERS
    element_validity.hpp
    lagrangeVector.hpp
    transMatrices.hpp
    cornerIndices.hpp
)

set(SOURCES
    lagrangeVector.cpp
    transMatrices.cpp
    cornerIndices.cpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(bezier PRIVATE ${SOURCES})
