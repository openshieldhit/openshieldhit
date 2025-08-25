# Create a library from local source files
# Cmake cannot track changes *.c glob files, so we list them explicitly
add_library(osh_particle
    osh_particle.c
    _osh_partdb.c
)

# Make sure consumers of the library see the headers
target_include_directories(osh_particle
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link math library
target_link_libraries(osh_particle m)
