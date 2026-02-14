cmake_minimum_required(VERSION 3.10)
project(openshieldhit C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# ---- Version Information ----
# Extract version from git at configure time and pass as compiler definition
find_package(Git QUIET)
if(GIT_FOUND)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --tags --always --dirty
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
endif()

# Fallback to project version if git is not available or not in a git repo
if(NOT GIT_VERSION)
    set(GIT_VERSION "0.0.1-unknown")
endif()

message(STATUS "OpenShieldHIT version: ${GIT_VERSION}")

# Compiler-specific warning flags
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-pedantic -Wall -Wextra)
elseif(MSVC)
    add_compile_options(/W4)
endif()

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
# Pass version as preprocessor definition (OSH_VERSION) to all targets
target_compile_definitions(osh_all INTERFACE
    OSH_VERSION="${GIT_VERSION}"
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
# ---- Installation and Packaging ----

# Install main executable
install(TARGETS osh
    RUNTIME DESTINATION bin
)

# Install example executables if built
if(OSH_BUILD_EXAMPLES)
    install(TARGETS gemca_sdl_viewer bnct_sdl
        RUNTIME DESTINATION bin
    )
    # Install example data files
    install(DIRECTORY examples/
        DESTINATION share/openshieldhit/examples
        FILES_MATCHING PATTERN "*.dat"
    )
endif()

# Install documentation
install(FILES README.md LICENSE
    DESTINATION share/openshieldhit
)

# ---- CPack Configuration ----
set(CPACK_PACKAGE_NAME "openshieldhit")
set(CPACK_PACKAGE_VENDOR "OpenShieldHIT Contributors")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Lean Monte Carlo Particle Transport in C")
set(CPACK_PACKAGE_DESCRIPTION "OpenShieldHIT is a modern Monte Carlo particle transport framework written entirely from scratch in C. It is designed to be lean, fast, efficient, readable, modular, and easy to extend.")
set(CPACK_PACKAGE_VERSION "${GIT_VERSION}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY ${CPACK_PACKAGE_NAME})

# Debian package specific settings
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "OpenShieldHIT Contributors")
set(CPACK_DEBIAN_PACKAGE_SECTION "science")
set(CPACK_DEBIAN_PACKAGE_HOMEPAGE "https://github.com/openshieldhit/openshieldhit")
# Optionally require libsdl2 if examples are built
if(OSH_BUILD_EXAMPLES)
    set(CPACK_DEBIAN_PACKAGE_DEPENDS "libsdl2-2.0-0 (>= 2.0.0)")
endif()

# RPM package specific settings (optional)
set(CPACK_RPM_PACKAGE_LICENSE "MIT")

# Generic tarball settings
set(CPACK_ARCHIVE_COMPONENT_INSTALL ON)

# Use DEB generator on Linux, others on respective platforms
if(UNIX AND NOT APPLE)
    set(CPACK_GENERATOR "DEB;TGZ")
elseif(APPLE)
    set(CPACK_GENERATOR "TGZ")
elseif(WIN32)
    set(CPACK_GENERATOR "ZIP")
endif()

include(CPack)