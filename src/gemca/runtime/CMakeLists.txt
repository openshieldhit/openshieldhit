add_library(osh_gemca_runtime
    osh_gemca_runtime.c
)

target_include_directories(osh_gemca_runtime
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

target_link_libraries(osh_gemca_runtime
    PRIVATE
        osh_gemca2
        osh_common
)

if(UNIX)
    target_link_libraries(osh_gemca_runtime PRIVATE m)
endif()
