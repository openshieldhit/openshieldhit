add_library(osh_neutron_physics
    osh_neutron_xsec.c
    osh_neutron_reaction.c
)

target_include_directories(osh_neutron_physics
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(osh_neutron_physics
    PRIVATE
        osh_common
        osh_nuclear
        $<$<BOOL:${UNIX}>:m>
)
