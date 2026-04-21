add_library(osh_dicom
    osh_dicom_parse.c
    osh_dicom_ct.c
    osh_dicom_rtdose.c
)

target_include_directories(osh_dicom
    PUBLIC
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(osh_dicom
    PRIVATE
        osh_common
)
