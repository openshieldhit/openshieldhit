add_library(osh_material
    osh_material.c
    osh_material_atomic_data.c
    osh_material_icru_db.c
    osh_material_loaddedx.c
    runtime/osh_material_prepare.c
)

target_include_directories(osh_material
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(osh_material
    PRIVATE
        osh_common
        osh_particle
        osh_physics
)
