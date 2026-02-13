cmake_minimum_required(VERSION 3.10)
project(openshieldhit C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

add_compile_options(-pedantic -Wall -Wextra)

# (Optional but recommended) default build type for Makefiles/Ninja
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING
        "Choose build type: Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif()

if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(
        $<$<CONFIG:Debug>:-Og>
        $<$<CONFIG:Debug>:-g3>
        $<$<CONFIG:Debug>:-fno-omit-frame-pointer>

        $<$<CONFIG:RelWithDebInfo>:-O3>
        $<$<CONFIG:RelWithDebInfo>:-g1>
        $<$<CONFIG:RelWithDebInfo>:-fno-omit-frame-pointer>

        $<$<CONFIG:Release>:-O3>
    )
endif()

add_compile_definitions($<$<CONFIG:Debug>:DEBUG>)

# ---- Libraries / modules (these create osh_common, osh_random, etc.) ----
add_subdirectory(src/common)
add_subdirectory(src/random)
add_subdirectory(src/beam)
add_subdirectory(src/particle)
add_subdirectory(src/gemca)
add_subdirectory(src/transport)

# ---- Umbrella target AFTER the libs exist ----
add_library(osh_all INTERFACE)
target_link_libraries(osh_all INTERFACE
    osh_common
    osh_random
    osh_beam
    osh_particle
    osh_gemca2
    osh_transport
)
target_include_directories(osh_all INTERFACE
    ${PROJECT_SOURCE_DIR}/src
)

# ---- Main executable ----
add_executable(osh src/main.c)
target_link_libraries(osh PRIVATE osh_all)

# ---- Tests ----
include(CTest)
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()

# ---- Examples ----
option(OSH_BUILD_EXAMPLES "Build example programs" ON)
if(OSH_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
