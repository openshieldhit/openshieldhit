add_library(osh_transport
    osh_transport.c
)

# Make sure consumers of the library see the headers
target_include_directories(osh_transport
    PUBLIC
    ${PROJECT_SOURCE_DIR}/src
)

# Link math library
target_link_libraries(osh_transport PRIVATE m)
