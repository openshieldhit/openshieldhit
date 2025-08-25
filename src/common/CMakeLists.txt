# Create a library from local source files
# Cmake cannot track changes *.c glob files, so we list them explicitly
add_library(osh_common
    osh_logger.c
    osh_rc.c
    osh_vect.c
    osh_file.c
    osh_readline.c
)

# Make sure consumers of the library see the headers
target_include_directories(osh_common
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link math library
target_link_libraries(osh_common m)
