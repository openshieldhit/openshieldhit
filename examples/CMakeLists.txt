find_package(SDL2 REQUIRED)

add_executable(gemca_sdl_viewer
    ${CMAKE_CURRENT_SOURCE_DIR}/gemca_sdl_viewer.c
)

target_link_libraries(gemca_sdl_viewer
    PRIVATE
        osh_gemca2
        osh_common
        osh_random
        osh_beam
        osh_particle
        SDL2::SDL2
        m
)


# Allow includes like "gemca/..." and "common/..."
target_include_directories(gemca_sdl_viewer
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)
