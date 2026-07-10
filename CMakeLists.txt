cmake_minimum_required(VERSION 3.10)
project(openshieldhit C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# ---- Version Information ----
# Extract version from git at configure time and pass as compiler definition.
# --match restricts the search to dotted semver release tags (vMAJOR.MINOR.PATCH,
# optionally with a suffix) so non-release tags (e.g. the test-fixtures-v1 data
# release, or a bare "v1") don't get picked up as the version and break the
# semver parsing below. The glob requires the two dots of a vX.Y.Z tag.
find_package(Git QUIET)
if(GIT_FOUND)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --tags --always --dirty --match "v[0-9]*.[0-9]*.[0-9]*"
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

# Parse semver components from tag (e.g. v1.2.3 or v1.2.3-5-gabcdef)
string(REGEX MATCH "([0-9]+)\\.([0-9]+)\\.([0-9]+)" OSH_SEMVER_MATCH "${GIT_VERSION}")
if(OSH_SEMVER_MATCH)
    set(OSH_VERSION_MAJOR "${CMAKE_MATCH_1}")
    set(OSH_VERSION_MINOR "${CMAKE_MATCH_2}")
    set(OSH_VERSION_PATCH "${CMAKE_MATCH_3}")
else()
    # Tagless checkout (shallow clone, fresh fork, CI mirror without tags):
    # `git describe --always` fell back to a bare commit hash with no X.Y.Z to
    # parse. Synthesise a semver-shaped default so OSH_VERSION always contains a
    # major.minor.patch that agrees with the component accessors, keeping the
    # commit as build metadata for provenance (e.g. "0.0.0+6aae2da").
    set(OSH_VERSION_MAJOR 0)
    set(OSH_VERSION_MINOR 0)
    set(OSH_VERSION_PATCH 0)
    set(GIT_VERSION "0.0.0+${GIT_VERSION}")
endif()
set(OSH_PROJECT_VERSION "${OSH_VERSION_MAJOR}.${OSH_VERSION_MINOR}.${OSH_VERSION_PATCH}")

message(STATUS "OpenShieldHIT version: ${GIT_VERSION}")

configure_file(
    ${PROJECT_SOURCE_DIR}/src/common/osh_version.h.in
    ${PROJECT_BINARY_DIR}/generated/common/osh_version.h
    @ONLY
)

# Compiler-specific warning flags
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-pedantic -Wall -Wextra)
elseif(MSVC)
    add_compile_options(/W4)
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS)   # fopen etc. are standard C; fopen_s is Annex K (optional)
    add_compile_definitions(_CRT_NONSTDC_NO_WARNINGS)  # strdup etc.: POSIX names kept for portability
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

set(OSH_VOXEL_LAYOUT "ROW_MAJOR" CACHE STRING
    "Default voxel storage layout compiled into the DCM/VOX path (ROW_MAJOR or MORTON8)")
set_property(CACHE OSH_VOXEL_LAYOUT PROPERTY STRINGS ROW_MAJOR MORTON8)
if(OSH_VOXEL_LAYOUT STREQUAL "ROW_MAJOR")
    add_compile_definitions(OSH_VOXEL_LAYOUT_DEFAULT=0)
elseif(OSH_VOXEL_LAYOUT STREQUAL "MORTON8")
    add_compile_definitions(OSH_VOXEL_LAYOUT_DEFAULT=8)
else()
    message(FATAL_ERROR "OSH_VOXEL_LAYOUT must be ROW_MAJOR or MORTON8")
endif()
message(STATUS "Voxel layout default: ${OSH_VOXEL_LAYOUT}")

option(OSH_ENABLE_COVERAGE "Enable GCC/Clang coverage instrumentation" OFF)
if(OSH_ENABLE_COVERAGE)
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
        add_compile_options(--coverage -O0 -g)
        if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
            add_compile_options(-fprofile-update=atomic)
        endif()
        add_link_options(--coverage)
    else()
        message(FATAL_ERROR "OSH_ENABLE_COVERAGE requires GCC or Clang")
    endif()
endif()

# Opt-in AddressSanitizer + UndefinedBehaviorSanitizer instrumentation. The
# `debug` preset is a plain -Og build; the sanitizer pass lives behind this
# switch and is exposed as the dedicated `asan` preset (build_asan/) so CI can
# run an instrumented ctest without slowing everyday debug builds. GCC/Clang get
# ASan+UBSan; MSVC only ships AddressSanitizer (no UBSan), so it gets ASan alone.
# -fno-sanitize-recover=all turns UBSan reports into hard aborts so a violation
# fails the test run instead of merely printing a warning.
option(OSH_ENABLE_SANITIZERS "Enable AddressSanitizer + UndefinedBehaviorSanitizer" OFF)
if(OSH_ENABLE_SANITIZERS)
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
        add_compile_options(-fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
        add_link_options(-fsanitize=address,undefined)
    elseif(MSVC)
        # MSVC ships AddressSanitizer only (no UBSan). Two MSVC-specific
        # requirements, both flagged in review: ASan is incompatible with the
        # /RTC run-time checks that CMake's Debug flags enable by default (strip
        # them, or cl.exe errors out), and with incremental linking (force it
        # off). The ASan runtime is import-linked automatically by cl.exe, but
        # its DLL (clang_rt.asan_dynamic-*.dll, shipped with the MSVC toolchain)
        # must be discoverable on PATH when the instrumented tests run.
        string(REGEX REPLACE "/RTC[1csu]+" "" CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}")
        add_compile_options(/fsanitize=address)
        add_link_options(/INCREMENTAL:NO)
    else()
        message(FATAL_ERROR "OSH_ENABLE_SANITIZERS requires GCC, Clang, or MSVC")
    endif()
endif()

# Opt-in native "quick-look" plot output (issue #238 prototype). OFF by default
# so the stock binary is byte-for-byte unchanged and carries zero plotting code;
# tools/ (matplotlib/pymchelper) stays the default plotting path. When enabled,
# the scoring library compiles a hand-rolled, dependency-free SVG line-plot
# writer reachable via `FileFormat SVG` on an Output block. The compile
# definition gates the dispatcher and the writer's translation unit alike.
option(OSH_ENABLE_PLOT "Enable native SVG plot output for scorer results" OFF)
if(OSH_ENABLE_PLOT)
    add_compile_definitions(OSH_ENABLE_PLOT)
    message(STATUS "Native SVG plot output: ENABLED (FileFormat SVG)")
endif()

# ---- Libraries / modules (these create osh_common, osh_random, etc.) ----
add_subdirectory(src/common)
add_subdirectory(src/random)
add_subdirectory(src/beam)
add_subdirectory(src/voxel)
add_subdirectory(src/material)
add_subdirectory(src/particle)
add_subdirectory(src/gemca)
add_subdirectory(src/scoring)
add_subdirectory(src/physics)
add_subdirectory(src/transport)
add_subdirectory(src/simulation)
add_subdirectory(src/dicom)
add_subdirectory(src/apps/osh)
add_subdirectory(src/apps/dicominspect)

# ---- Umbrella target AFTER the libs exist ----
add_library(osh_all INTERFACE)
target_link_libraries(osh_all INTERFACE
    osh_common
    osh_random
    osh_beam
    osh_beam_runtime
    osh_material
    osh_particle
    osh_gemca2
    osh_gemca_runtime
    osh_scoring
    osh_transport
    osh_dicom
)
target_include_directories(osh_all INTERFACE
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_BINARY_DIR}/generated
)

# ---- Internal CLI parser (used by executable + tests, not installed) ----
add_library(osh_cli STATIC src/cli/osh_cli.c)
target_include_directories(osh_cli PUBLIC ${PROJECT_SOURCE_DIR}/src)
# osh_common provides the shared duration parser used by --max-time.
target_link_libraries(osh_cli PUBLIC osh_common)

# ---- Main executable ----
add_executable(openshieldhit src/apps/osh/main.c)
target_link_libraries(openshieldhit PRIVATE osh_app_osh osh_all osh_cli)

# ---- Tests ----
include(CTest)
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()

# ---- Examples ----
option(OSH_BUILD_EXAMPLES "Build example programs" OFF)
if(OSH_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
# ---- Installation and Packaging ----

# Install main executable
install(TARGETS openshieldhit
    RUNTIME DESTINATION bin
)

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
# Sanitize version for Debian (must start with digit per Debian policy)
# Strip leading 'v' and prepend '0.' if version starts with non-digit
string(REGEX REPLACE "^v" "" DEBIAN_VERSION "${GIT_VERSION}")
if(NOT DEBIAN_VERSION MATCHES "^[0-9]")
    # Version starts with non-digit (e.g., commit hash), prepend "0." for valid Debian version
    set(DEBIAN_VERSION "0.${DEBIAN_VERSION}")
endif()
set(CPACK_DEBIAN_PACKAGE_VERSION "${DEBIAN_VERSION}")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "OpenShieldHIT Contributors")
set(CPACK_DEBIAN_PACKAGE_SECTION "science")
set(CPACK_DEBIAN_PACKAGE_HOMEPAGE "https://github.com/openshieldhit/openshieldhit")

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
