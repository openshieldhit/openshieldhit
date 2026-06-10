add_library(osh_nuclear
    osh_nuclear_tripathi.c
    osh_nuclear_pp.c
    osh_nuclear_abrasion.c
    osh_nuclear_handler.c
)

target_include_directories(osh_nuclear
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(osh_nuclear
    PRIVATE
        osh_kinematics
        osh_random
        $<$<BOOL:${UNIX}>:m>
)
