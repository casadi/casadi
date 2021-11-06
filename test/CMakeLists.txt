add_subdirectory(ref)
add_subdirectory(drivers)

# cmake_policy(SET CMP0094 NEW)
# set(Python3_FIND_STRATEGY "LOCATION")
find_package(Python3 REQUIRED COMPONENTS Interpreter)
message(STATUS "Found Python ${Python3_VERSION}")

set(PYTHON_TEST_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/python)
set(PYTHON_TEST_GEN_FILE
    ${PYTHON_TEST_DIR}/gen-quadratic-optimal-control-problem-test.py)
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/gen/test_alm_mcp_24x24.cpp
    COMMAND Python3::Interpreter 
        ${PYTHON_TEST_GEN_FILE}
        # --nx 40 --nu 40 --seed 55
        --nx 24 --nu 24 --seed 55
        --output ${CMAKE_CURRENT_BINARY_DIR}/gen --name test_alm_mcp_24x24
    MAIN_DEPENDENCY ${PYTHON_TEST_GEN_FILE}
)
# set(gen_test_files ${CMAKE_CURRENT_BINARY_DIR}/gen/test_alm_mcp_24x24.cpp)
set(gen_test_files "")

# Test executable compilation and linking
file(GLOB test_files *.cpp)
add_executable(tests ${test_files} ${gen_test_files})
target_include_directories(tests 
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}
)
target_link_libraries(tests
    PRIVATE
        GTest::gtest_main
        GTest::gmock
        alpaqa::alpaqa
        alpaqa-ref
)

# Add tests
if (NOT CMAKE_CROSSCOMPILING)
    gtest_discover_tests(tests)
endif()