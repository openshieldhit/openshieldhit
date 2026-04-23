add_library(osh_gemca_runtime
    osh_gemca_runtime.c
    osh_gemca_runtime_voxel.c
)

target_include_directories(osh_gemca_runtime
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(osh_gemca_runtime
    PRIVATE
        osh_gemca2
        osh_common
        osh_raytrace
)

if(UNIX)
    target_link_libraries(osh_gemca_runtime PRIVATE m)
endif()

# --- Optional AVX2 + FMA acceleration ---
include(CheckCCompilerFlag)
check_c_compiler_flag("-mavx2" GEMCA_HAVE_AVX2)
check_c_compiler_flag("-mfma"  GEMCA_HAVE_FMA)

if(GEMCA_HAVE_AVX2 AND GEMCA_HAVE_FMA)
    target_sources(osh_gemca_runtime PRIVATE osh_gemca_runtime_avx2.c)
    set_source_files_properties(osh_gemca_runtime_avx2.c
        PROPERTIES COMPILE_OPTIONS "-mavx2;-mfma"
    )
    target_compile_definitions(osh_gemca_runtime PRIVATE OSH_GEMCA_RUNTIME_HAVE_AVX2=1)
    message(STATUS "gemca_runtime: AVX2+FMA zone batch enabled")
else()
    message(STATUS "gemca_runtime: AVX2+FMA not available, using scalar fallback")
endif()
