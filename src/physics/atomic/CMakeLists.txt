add_library(osh_physics
    osh_physics_bethe.c
    osh_physics_scat.c
    osh_physics_scat_highland.c
    osh_physics_scat_moliere.c
    osh_physics_straggling.c
)

target_include_directories(osh_physics
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(osh_physics
    PRIVATE
        osh_common
        osh_kinematics
        $<$<BOOL:${UNIX}>:m>
)
