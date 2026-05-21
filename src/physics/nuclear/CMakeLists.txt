add_library(osh_nuclear
    osh_nuclear_tripathi.c
)

target_include_directories(osh_nuclear
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(osh_nuclear
    PRIVATE
        $<$<BOOL:${UNIX}>:m>
)
