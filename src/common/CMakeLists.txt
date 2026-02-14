# Create a library from local source files
# Cmake cannot track changes *.c glob files, so we list them explicitly
add_library(osh_common
    osh_logger.c
    osh_rc.c
    osh_vect.c
    osh_file.c
    osh_readline.c
    osh_coord.c
)

# Make sure consumers of the library see the headers
target_include_directories(osh_common
    PUBLIC
    ${PROJECT_SOURCE_DIR}/src
)

# Link math library (Unix/Linux only - Windows has it in the C runtime)
if(UNIX)
    target_link_libraries(osh_common PRIVATE m)
endif()
