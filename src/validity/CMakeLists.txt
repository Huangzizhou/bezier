set(HEADERS
	ValidityChecker.hpp
	cornerIndices.hpp
	lagrangeVector.hpp
	transMatrices.hpp
)

set(SOURCES
	ValidityChecker.cpp
	cornerIndices.cpp
)

add_subdirectory(lagvec)
add_subdirectory(lagvecT)
add_subdirectory(transmat)
add_subdirectory(transmatT)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(bezier PRIVATE ${SOURCES})