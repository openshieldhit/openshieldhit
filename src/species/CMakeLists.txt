add_library(osh_species
    osh_isotope_db_generated.c
)

# Make sure consumers of the library see the headers
target_include_directories(osh_species
    PUBLIC
    ${PROJECT_SOURCE_DIR}/src
)
