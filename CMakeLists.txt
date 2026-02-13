cmake_minimum_required(VERSION 3.10)
project(openshieldhit C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

add_compile_options(-pedantic -Wall -Wextra)

if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(
        $<$<CONFIG:Debug>:-Og>
        $<$<CONFIG:Debug>:-g1>
        $<$<CONFIG:Debug>:-fno-omit-frame-pointer>

        $<$<CONFIG:Release>:-O3>
    )
endif()
# Enable DEBUG macro in Debug builds
add_compile_definitions($<$<CONFIG:Debug>:DEBUG>)

# Libraries / modules
add_subdirectory(src/common)
add_subdirectory(src/random)
add_subdirectory(src/beam)
add_subdirectory(src/particle)
add_subdirectory(src/gemca)
add_subdirectory(src/transport)


# Main executable
add_executable(osh src/main.c)

# Link internal libraries
target_link_libraries(osh
    PRIVATE
        osh_common
        osh_beam
        osh_particle
        osh_random
        osh_gemca2
        osh_transport
)


# Tests (canonical)
include(CTest)          # defines BUILD_TESTING + calls enable_testing()
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()

# Options
option(OSH_BUILD_EXAMPLES "Build example programs" ON)

# Examples (optional)
if(OSH_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()