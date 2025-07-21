cmake_minimum_required(VERSION 3.10)
project(openshieldhit C)

set(CMAKE_C_STANDARD 90)
set(CMAKE_C_STANDARD_REQUIRED ON)
add_compile_options(-pedantic -Wall -Wextra)

# Global include directory
include_directories(${CMAKE_SOURCE_DIR}/include)

# Add subdirectories
add_subdirectory(src/common)

enable_testing()
add_subdirectory(tests)

# Main executable
add_executable(osh src/main.c)
target_link_libraries(osh PRIVATE osh_common)


