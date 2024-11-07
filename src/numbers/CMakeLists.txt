set(SOURCES
    Interval.cpp
    Matrix.cpp
    Rational.cpp
)

set(HEADERS
    Interval.hpp
    Rational.hpp
    Matrix.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(bezier PRIVATE ${SOURCES})
