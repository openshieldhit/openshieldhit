add_library(osh_random
    osh_rng_pcg32.c
    osh_rng_xoshiro256ss.c
    osh_rng.c
)

# Make sure consumers of the library see the headers
target_include_directories(osh_random
    PUBLIC
    ${PROJECT_SOURCE_DIR}/src
)

# Link math library
target_link_libraries(osh_random PRIVATE m)
