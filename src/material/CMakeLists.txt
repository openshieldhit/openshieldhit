add_library(osh_material
    osh_material.c
    osh_material_parse.c
)

target_include_directories(osh_material
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(osh_material
    PRIVATE
        osh_common
)
