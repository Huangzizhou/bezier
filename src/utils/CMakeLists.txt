set(HEADERS
    globals.hpp
    Timer.hpp
    combinatorics.hpp
    Settings.hpp
    helpString.hpp
)
set(SOURCES
    Settings.cpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(bezier PRIVATE ${SOURCES})
