cmake_minimum_required(VERSION 3.10)
project(wavecar_compressor LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE ON)

# Add all .cpp files in src/
file(GLOB SOURCES CONFIGURE_DEPENDS src/*.cpp)

add_executable(wavecar_compressor ${SOURCES})

target_include_directories(wavecar_compressor
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src
)

add_custom_target(
    sca
    COMMAND clang-tidy -p ${CMAKE_BINARY_DIR} ${SOURCES} -checks=* -warnings-as-errors=*
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Running static code analysis"
)


# Optional: add warnings
target_compile_options(wavecar_compressor PRIVATE -Wall -Wextra -Wpedantic)
