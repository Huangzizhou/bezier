set(SOURCES
	lagrangeVectorT_1_1_1.cpp
	lagrangeVectorT_1_1_2.cpp
	lagrangeVectorT_1_1_3.cpp
	lagrangeVectorT_1_1_4.cpp
	lagrangeVectorT_1_1_5.cpp
	lagrangeVectorT_2_1_1.cpp
	lagrangeVectorT_2_1_2.cpp
	lagrangeVectorT_2_2_1.cpp
	lagrangeVectorT_2_2_2.cpp
	lagrangeVectorT_2_2_3.cpp
	lagrangeVectorT_3_1_1.cpp
	lagrangeVectorT_3_3_1.cpp
	lagrangeVectorT_3_3_2.cpp
	lagrangeVectorT_2_2_4.cpp
	lagrangeVectorT_3_3_3.cpp
	# lagrangeVectorT_3_3_4.cpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(bezier PRIVATE ${SOURCES})