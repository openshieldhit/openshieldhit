set(OSH_RAYTRACE_ALGORITHM "SIDDON" CACHE STRING
    "Voxel raytrace algorithm. SIMPLE = 3-D DDA pedagogical reference; SIDDON = Siddon parametric with bounding-box clip (default); JACOBS = alpha-parametric Jacobs (1998).")
set_property(CACHE OSH_RAYTRACE_ALGORITHM PROPERTY STRINGS SIMPLE SIDDON JACOBS)

add_library(osh_raytrace STATIC)

if(OSH_RAYTRACE_ALGORITHM STREQUAL "SIMPLE")
    target_sources(osh_raytrace PRIVATE osh_raytrace_simple_msh.c)
elseif(OSH_RAYTRACE_ALGORITHM STREQUAL "SIDDON")
    target_sources(osh_raytrace PRIVATE osh_raytrace_siddon_msh.c)
elseif(OSH_RAYTRACE_ALGORITHM STREQUAL "JACOBS")
    target_sources(osh_raytrace PRIVATE osh_raytrace_jacobs_msh.c)
else()
    message(FATAL_ERROR
        "Unknown OSH_RAYTRACE_ALGORITHM: '${OSH_RAYTRACE_ALGORITHM}'. "
        "Valid options: SIMPLE, SIDDON, JACOBS.")
endif()

# Expose the active algorithm as a compile definition so diagnostic code
# can report which implementation is in use.
target_compile_definitions(osh_raytrace
    PUBLIC OSH_RAYTRACE_${OSH_RAYTRACE_ALGORITHM}=1)

target_include_directories(osh_raytrace
    PUBLIC ${PROJECT_SOURCE_DIR}/src)

target_link_libraries(osh_raytrace PRIVATE $<$<BOOL:${UNIX}>:m>)
