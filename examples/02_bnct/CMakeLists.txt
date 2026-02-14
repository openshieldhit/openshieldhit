find_package(SDL2 REQUIRED)

add_executable(bnct_sdl
        ${CMAKE_CURRENT_SOURCE_DIR}/bnct_sdl.c
)


target_link_libraries(bnct_sdl
    PRIVATE
        osh_all
        SDL2::SDL2
)

# Link math library (Unix/Linux only - Windows has it in the C runtime)
if(UNIX)
    target_link_libraries(bnct_sdl PRIVATE m)
endif()
