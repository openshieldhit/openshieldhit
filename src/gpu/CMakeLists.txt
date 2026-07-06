# src/gpu/ — experimental CUDA backend (built only when OSH_ENABLE_CUDA=ON).
#
# Layering (DEVELOPER.md §9): this module consumes the runtime layers; no
# CPU-path target may depend on it.  Kernels contain no physics of their
# own — they call the OSH_HD static-inline bodies shared with the CPU.

add_library(osh_gpu STATIC
    osh_gpu_gemca_mirror.cu
)

target_include_directories(osh_gpu
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_BINARY_DIR}/generated
)

target_link_libraries(osh_gpu
    PUBLIC
        osh_gemca_runtime
        CUDA::cudart
)

set_target_properties(osh_gpu PROPERTIES
    CUDA_ARCHITECTURES "80"          # A100 (Athena); extend to "80;90" for Hopper/GH200
    CUDA_SEPARABLE_COMPILATION OFF
)

# Device code gets its warning story here (the top-level C flags are scoped
# to $<COMPILE_LANGUAGE:C> and never reach nvcc).
target_compile_options(osh_gpu PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:--Werror=all-warnings>
)

# CPU/GPU zone-membership parity + throughput harness.
add_executable(osh_gpu_zone_parity osh_gpu_zone_parity_main.c)
target_link_libraries(osh_gpu_zone_parity PRIVATE
    osh_gpu
    osh_app_osh
    osh_gemca2
    osh_common
)

# RNG device known-answer + throughput; also the compile guard that keeps
# the RNG _hd headers genuinely device-compilable.
add_executable(osh_gpu_rng_bench osh_gpu_rng_bench.cu)
target_include_directories(osh_gpu_rng_bench PRIVATE ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(osh_gpu_rng_bench PRIVATE osh_random)
set_target_properties(osh_gpu_rng_bench PROPERTIES CUDA_ARCHITECTURES "80")

# FP64 atomicAdd contention on scoring-like bin patterns (issue #231 D4).
add_executable(osh_gpu_atomics_bench osh_gpu_atomics_bench.cu)
target_include_directories(osh_gpu_atomics_bench PRIVATE ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/include)
set_target_properties(osh_gpu_atomics_bench PROPERTIES CUDA_ARCHITECTURES "80")

# Device-buildability guard: compiles the entire OSH_HD surface as device
# code (rule R2, docs/dev/gpu_port_plan.md).
add_executable(osh_gpu_hd_guard osh_gpu_hd_guard.cu)
target_include_directories(osh_gpu_hd_guard PRIVATE ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/include)
set_target_properties(osh_gpu_hd_guard PROPERTIES CUDA_ARCHITECTURES "80")
