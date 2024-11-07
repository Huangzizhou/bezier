set(SOURCES
    Interval.cpp
    Matrix.cpp
)

set(HEADERS
    Interval.hpp
    Rational.hpp
    Matrix.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(bezier PRIVATE ${SOURCES})
