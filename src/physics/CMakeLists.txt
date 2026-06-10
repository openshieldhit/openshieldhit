add_library(osh_kinematics
    osh_kinematics.c
)

target_include_directories(osh_kinematics
    PUBLIC
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(osh_kinematics
    PRIVATE
        osh_common
        osh_random
        $<$<BOOL:${UNIX}>:m>
)

add_subdirectory(atomic)
add_subdirectory(nuclear)
