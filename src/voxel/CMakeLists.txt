add_library(osh_voxel
    osh_voxel_hu_lut.c
    osh_voxel_hu_lut_permatassari.c
)

# Allow includes like "voxel/osh_voxel_hu_lut.h" (src/) and
# "openshieldhit/voxel.h" (include/).
target_include_directories(osh_voxel
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(osh_voxel
    PRIVATE
        osh_common
)
