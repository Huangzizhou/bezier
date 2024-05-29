set(SOURCES
    # IndexPermutator.cpp
)

set(HEADERS
    # IndexPermutator.hpp
    # SubdivisionHierarchy.hpp
    Timer.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(bezier PUBLIC ${SOURCES})
