set(HEADERS
	element_validity.hpp
	cornerIndices.hpp
	lagrangeVector.hpp
	transMatrices.hpp
)

set(SOURCES
	transMatrices_1_1_1.cpp
	lagrangeVector_1_1_1.cpp
	transMatrices_1_1_2.cpp
	lagrangeVector_1_1_2.cpp
	transMatrices_1_1_3.cpp
	lagrangeVector_1_1_3.cpp
	transMatrices_1_1_4.cpp
	lagrangeVector_1_1_4.cpp
	transMatrices_1_1_5.cpp
	lagrangeVector_1_1_5.cpp
	transMatrices_2_1_1.cpp
	lagrangeVector_2_1_1.cpp
	transMatrices_2_1_2.cpp
	lagrangeVector_2_1_2.cpp
	transMatrices_2_2_1.cpp
	lagrangeVector_2_2_1.cpp
	transMatrices_2_2_2.cpp
	lagrangeVector_2_2_2.cpp
	transMatrices_2_2_3.cpp
	lagrangeVector_2_2_3.cpp
	transMatrices_2_2_4.cpp
	lagrangeVector_2_2_4.cpp
	transMatrices_3_1_1.cpp
	lagrangeVector_3_1_1.cpp
	transMatrices_3_1_2.cpp
	lagrangeVector_3_1_2.cpp
	transMatrices_3_3_1.cpp
	lagrangeVector_3_3_1.cpp
	transMatrices_3_3_2.cpp
	lagrangeVector_3_3_2.cpp
	transMatrices_3_3_3.cpp
	lagrangeVector_3_3_3.cpp
	transMatrices_3_3_4.cpp
	lagrangeVector_3_3_4.cpp
	cornerIndices.cpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(bezier PRIVATE ${SOURCES})