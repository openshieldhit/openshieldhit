add_library(osh_material
    osh_material.c
    osh_material_atomic_data.c
    osh_material_icru_db.c
    osh_material_loaddedx.c
    osh_material_parse.c
)

target_include_directories(osh_material
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(osh_material
    PRIVATE
        osh_common
        osh_particle
)
