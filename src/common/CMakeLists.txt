# Create a library from local source files
add_library(osh_common
    osh_logger.c
    osh_vect.c
)

# Make sure consumers of the library see the headers
target_include_directories(osh_common
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)

# Link math library
target_link_libraries(osh_common m)
