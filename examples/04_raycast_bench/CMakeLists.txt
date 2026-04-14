find_package(SDL2 REQUIRED)

add_executable(gemca_raycast_bench
    ${CMAKE_CURRENT_SOURCE_DIR}/gemca_raycast_bench.c
)


target_link_libraries(gemca_raycast_bench
    PRIVATE
        osh_all
        SDL2::SDL2
)

if(UNIX)
    target_link_libraries(gemca_raycast_bench PRIVATE m)
endif()
