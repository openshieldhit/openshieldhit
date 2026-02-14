find_package(SDL2 REQUIRED)

add_executable(gemca_sdl_viewer
    ${CMAKE_CURRENT_SOURCE_DIR}/gemca_sdl_viewer.c
)

target_link_libraries(gemca_sdl_viewer
    PRIVATE
        osh_all
        SDL2::SDL2
)

# Link math library (Unix/Linux only - Windows has it in the C runtime)
if(UNIX)
    target_link_libraries(gemca_sdl_viewer PRIVATE m)
endif()
