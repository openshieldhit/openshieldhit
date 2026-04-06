add_executable(gemca_bench
    ${CMAKE_CURRENT_SOURCE_DIR}/gemca_bench.c
)

target_link_libraries(gemca_bench
    PRIVATE
        osh_all
)

if(UNIX)
    target_link_libraries(gemca_bench PRIVATE m)
endif()
