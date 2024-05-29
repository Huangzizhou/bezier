set(SOURCES
    globals.hpp
    # element_validity.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(jacobian PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################

add_subdirectory(utils)
add_subdirectory(numbers)
# add_subdirectory(algorithms)
