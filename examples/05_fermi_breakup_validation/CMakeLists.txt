add_executable(fbu_scan
    ${CMAKE_CURRENT_SOURCE_DIR}/fbu_scan.c
)

target_include_directories(fbu_scan PRIVATE
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(fbu_scan
    PRIVATE
        osh_nuclear
        osh_particle
        osh_kinematics
        osh_random
        osh_common
)

if(UNIX)
    target_link_libraries(fbu_scan PRIVATE m)
endif()
