cmake_minimum_required(VERSION 3.10)

project(smolpack
    VERSION 0.1.0
    LANGUAGES C
)

# ----------------------------
# Library sources
# ----------------------------

set(SMOLPACK_SRC
    src/smolpack.c
    src/smolyak_init.c
    src/smolyak_weights.c
    src/smolyak_alg.c
    src/cc_init.c
    src/cc_weights.c
    src/cc_alg.c
)

# ----------------------------
# Static library
# ----------------------------

add_library(smolpack STATIC ${SMOLPACK_SRC})

target_include_directories(smolpack
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:include>
)

set_target_properties(smolpack PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 0
)

# ----------------------------
# Shared library
# ----------------------------

add_library(smolpack_shared SHARED ${SMOLPACK_SRC})

target_include_directories(smolpack_shared
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:include>
)

set_target_properties(smolpack_shared PROPERTIES
    OUTPUT_NAME smolpack
    VERSION ${PROJECT_VERSION}
    SOVERSION 0
)

# ----------------------------
# Installation
# ----------------------------

install(
    FILES src/smolpack.h
    DESTINATION include
)

install(
    TARGETS smolpack smolpack_shared
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)

# ----------------------------
# Tests
# ----------------------------

option(SMOLPACK_BUILD_TESTS "Build C tests" ON)

if(SMOLPACK_BUILD_TESTS)

    enable_testing()

    add_executable(test_smolpack
        tests/test_smolpack.c
        tests/genz.c
    )

    target_link_libraries(test_smolpack
        smolpack
        $<$<NOT:$<PLATFORM_ID:Windows>>:m>
    )

    add_test(
        NAME test_smolpack
        COMMAND ${CMAKE_BINARY_DIR}/test_smolpack
    )

endif()

message(STATUS "SMOLPACK C library CMake setup complete.")