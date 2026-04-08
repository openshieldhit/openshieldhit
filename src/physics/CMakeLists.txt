add_library(osh_physics
    osh_physics_bethe.c
)

target_include_directories(osh_physics
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(osh_physics
    PRIVATE
        osh_common
        $<$<BOOL:${UNIX}>:m>
)
