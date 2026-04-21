add_executable(dicominspect main.c)

target_include_directories(dicominspect PRIVATE
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(dicominspect PRIVATE
    osh_common
    osh_dicom
)
