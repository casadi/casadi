# Test executable compilation and linking
file(GLOB test_files *.cpp)
add_executable(tests ${test_files})
target_include_directories(tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(tests
    PRIVATE
        GTest::gtest_main
        GTest::gmock
        panoc-alm::panoc-alm
)

# Add tests
if(NOT CMAKE_CROSSCOMPILING)
    gtest_discover_tests(tests)
endif()