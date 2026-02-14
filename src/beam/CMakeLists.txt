# Create a library from local source files
# Cmake cannot track changes *.c glob files, so we list them explicitly
add_library(osh_beam
    osh_beam.c
    osh_beam_spots.c
)

# Make sure consumers of the library see the headers
target_include_directories(osh_beam
    PUBLIC
    ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(osh_beam
    PRIVATE
        osh_common
        osh_particle
)

# Link math library (Unix/Linux only - Windows has it in the C runtime)
if(UNIX)
    target_link_libraries(osh_beam PRIVATE m)
endif()