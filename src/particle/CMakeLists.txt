add_library(osh_particle
    osh_particle.c
    osh_particle_db.c
    osh_isotope_db.c
    osh_isotope_db_generated.c
)

# Make sure consumers of the library see the headers
target_include_directories(osh_particle
    PUBLIC
    ${PROJECT_SOURCE_DIR}/src
)
