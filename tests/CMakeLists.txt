add_executable(test_osh_vect
    test_osh_vect.c
)

target_include_directories(test_osh_vect
    PRIVATE ${PROJECT_SOURCE_DIR}/include
)

# Link to the common code
target_link_libraries(test_osh_vect
    PRIVATE osh_common
)

add_test(NAME osh_vect_tests COMMAND test_osh_vect)

# --- Add logger test below ---



add_executable(test_osh_logger
    test_osh_logger.c
)

target_include_directories(test_osh_logger
    PRIVATE ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(test_osh_logger
    PRIVATE osh_common
)

add_test(NAME osh_logger_tests COMMAND test_osh_logger)