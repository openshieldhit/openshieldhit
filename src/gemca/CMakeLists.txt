add_subdirectory(runtime)

# Gemca2 library
add_library(osh_gemca2
    osh_gemca2.c
    osh_gemca2_calc_body.c
    osh_gemca2_calc_surface.c
    osh_gemca2_calc_zone.c
    osh_gemca2_dist.c
    osh_geometry_workspace.c

    prepare/osh_gemca_body_prepare.c
    prepare/osh_gemca_zone_compile.c
    prepare/osh_gemca_stack.c

    voxel/osh_gemca2_voxel_hu.c
    voxel/osh_gemca2_voxel_hu_permatassari.c
)

# Allow includes like "gemca/..." and "common/..." (src/) and
# "openshieldhit/geometry.h" (include/).
target_include_directories(osh_gemca2
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/include
)

# common headers/functions, link it
target_link_libraries(osh_gemca2
    PRIVATE
        osh_common
        osh_material
)

# Link math library (Unix/Linux only - Windows has it in the C runtime)
if(UNIX)
    target_link_libraries(osh_gemca2 PRIVATE m)
endif()
