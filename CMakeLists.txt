cmake_minimum_required(VERSION 3.10)
project(openshieldhit C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
add_compile_options(-pedantic -Wall -Wextra)

# Enable DEBUG macro in Debug builds
add_compile_definitions($<$<CONFIG:Debug>:DEBUG>)

# Add subdirectories
add_subdirectory(src/common)
add_subdirectory(src/beam)
add_subdirectory(src/particle)

enable_testing()
add_subdirectory(tests)

# Main executable
add_executable(osh src/main.c)

# Link internal libraries
target_link_libraries(osh
    PRIVATE
        osh_common
        osh_beam
        osh_particle
)


