set(OSH_RAYTRACE_ALGORITHM "SIDDON" CACHE STRING
    "Voxel raytrace algorithm. SIMPLE = 3-D DDA pedagogical reference; SIDDON = Siddon parametric with bounding-box clip (default).")
set_property(CACHE OSH_RAYTRACE_ALGORITHM PROPERTY STRINGS SIMPLE SIDDON)

add_library(osh_raytrace STATIC)

if(OSH_RAYTRACE_ALGORITHM STREQUAL "SIMPLE")
    target_sources(osh_raytrace PRIVATE osh_raytrace_simple.c)
elseif(OSH_RAYTRACE_ALGORITHM STREQUAL "SIDDON")
    target_sources(osh_raytrace PRIVATE osh_raytrace_siddon.c)
else()
    message(FATAL_ERROR
        "Unknown OSH_RAYTRACE_ALGORITHM: '${OSH_RAYTRACE_ALGORITHM}'. "
        "Valid options: SIMPLE, SIDDON.")
endif()

# Expose the active algorithm as a compile definition so diagnostic code
# can report which implementation is in use.
target_compile_definitions(osh_raytrace
    PUBLIC OSH_RAYTRACE_${OSH_RAYTRACE_ALGORITHM}=1)

target_include_directories(osh_raytrace
    PUBLIC ${PROJECT_SOURCE_DIR}/src)

target_link_libraries(osh_raytrace PRIVATE $<$<BOOL:${UNIX}>:m>)
