# Vendored MCPL core (see README.md).  mcpl.c includes <zlib.h> unconditionally, so
# this target needs zlib to link even though osh_scoring_save_mcpl.c never calls a
# gz* function itself.  Prefer a system zlib; fall back to fetching and building
# madler/zlib from source so a fresh checkout builds without any manual dependency
# install, on any of the three CI platforms.
find_package(ZLIB QUIET)
if(NOT ZLIB_FOUND)
    message(STATUS "System zlib not found; fetching madler/zlib for the vendored MCPL core")
    include(FetchContent)
    set(SKIP_INSTALL_ALL ON CACHE BOOL "" FORCE)
    FetchContent_Declare(
        osh_zlib
        GIT_REPOSITORY https://github.com/madler/zlib.git
        GIT_TAG v1.3.1
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(osh_zlib)
    # Upstream zlib's CMakeLists does not mark zconf.h's generated location (the
    # build tree) as a usage requirement of zlibstatic, so a consumer outside the
    # zlib subdirectory scope cannot find it without this. Same for the source
    # tree's own zlib.h.
    target_include_directories(zlibstatic PUBLIC
        $<BUILD_INTERFACE:${osh_zlib_SOURCE_DIR}>
        $<BUILD_INTERFACE:${osh_zlib_BINARY_DIR}>
    )
    add_library(ZLIB::ZLIB ALIAS zlibstatic)
endif()

add_library(osh_thirdparty_mcpl STATIC
    mcpl.c
    mcpl_fileutils.c
)
target_include_directories(osh_thirdparty_mcpl PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(osh_thirdparty_mcpl PRIVATE ZLIB::ZLIB)

# mcpl.h's MCPL_API expands to __declspec(dllimport) on Windows unless this is
# defined, since upstream defaults to being consumed as a DLL. We always build
# a plain static library (see the README.md rationale), so every mcpl.h
# include — this target's own .c files and any consumer's, e.g.
# osh_scoring_save_mcpl.c — needs the same define; PUBLIC propagates it to
# whatever links this target (privately or not), matching how the include
# directory above is already propagated.
target_compile_definitions(osh_thirdparty_mcpl PUBLIC mcpl_LIB_IS_STATIC)

# Vendored third-party source: not held to this project's own warning/format rules.
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(osh_thirdparty_mcpl PRIVATE -w)
elseif(MSVC)
    target_compile_options(osh_thirdparty_mcpl PRIVATE /W0)
endif()

# mcpl.c serialises a particle record by casting a byte offset into a
# char[MCPLIMP_MAX_PARTICLE_SIZE] buffer to double*/float*/int32_t* — a known,
# upstream-accepted pattern that UBSan's alignment check flags (the buffer sits
# inside a heap-allocated struct at a byte offset the compiler is free to leave
# non-8-aligned, since a char[] member itself needs only 1-byte alignment).
# Harmless on every platform this project targets (x86-64, ARM64 both tolerate
# misaligned double/int32_t loads and stores), so scope just the alignment
# check off for this one vendored target rather than patch upstream source or
# silence -fsanitize=undefined project-wide.
if(OSH_ENABLE_SANITIZERS AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(osh_thirdparty_mcpl PRIVATE -fno-sanitize=alignment)
endif()
