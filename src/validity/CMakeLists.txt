set(HEADERS
	Validator.hpp
	StaticValidator.hpp
	ContinuousValidator.hpp
	cornerIndices.hpp
	lagrangeVector.hpp
	lagrangeEvaluate.hpp
	transMatrices.hpp
)

set(SOURCES
	Validator.cpp
	cornerIndices.cpp
)

add_subdirectory(lagvec)
add_subdirectory(lagvecT)
add_subdirectory(transmat)
add_subdirectory(transmatT)
add_subdirectory(lageval)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(bezier PRIVATE ${SOURCES})