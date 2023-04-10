cmake_minimum_required(VERSION 3.21)
project(alpaqa-test C CXX)
enable_testing()
include(GoogleTest)

# Dependencies
find_package(GTest REQUIRED)
if (NOT TARGET alpaqa::alpaqa)
    find_package(alpaqa REQUIRED)
endif()

add_executable(tests
    "accelerators/test-lbfgs.cpp"
    "inner/test-panoc.cpp"
    "util/test-type-erasure.cpp"
    "util/test-index-set.cpp"
    "util/test-set-intersection.cpp"
    "util/test-sparse-ops.cpp"
    "util/io/test-csv.cpp"
    "outer/test-alm.cpp"
    "problem/test-type-erased-problem.cpp"
)
target_include_directories(tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(tests PRIVATE _CRT_SECURE_NO_WARNINGS)
target_link_libraries(tests PRIVATE alpaqa::alpaqa alpaqa::warnings
                                    GTest::gtest_main 
                                    GTest::gmock)

if (TARGET alpaqa::casadi-loader)
    target_sources(tests PUBLIC "outer/test-casadi.cpp")
    add_library(rosenbrock_functions_test MODULE
        "outer/rosenbrock_functions_test.c")
    set_target_properties("rosenbrock_functions_test" PROPERTIES
        PREFIX ""
        RELEASE_POSTFIX ""
        DEBUG_POSTFIX ""
        RELWITHDEBINFO_POSTFIX ""
        MINSIZEREL_POSTFIX ""
        OUTPUT_NAME "rosenbrock_functions_test")
    target_compile_definitions(tests PRIVATE
        ROSENBROCK_FUNC_DLL=\"$<TARGET_FILE:rosenbrock_functions_test>\")
    add_dependencies(tests rosenbrock_functions_test)
    target_link_libraries(tests PRIVATE alpaqa::casadi-loader)
    add_custom_command(TARGET rosenbrock_functions_test POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/outer/rosenbrock_functions_test.csv
                $<TARGET_FILE_DIR:rosenbrock_functions_test>)
endif()

if (NOT CMAKE_CROSSCOMPILING)
    gtest_discover_tests(tests DISCOVERY_TIMEOUT 60)
endif()
add_executable(alpaqa::tests ALIAS tests)

if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.21) 
    # Copy DLLs for test executable
    add_custom_command(TARGET tests PRE_LINK COMMAND
        COMMAND ${CMAKE_COMMAND} -E 
            $<IF:$<BOOL:$<TARGET_RUNTIME_DLLS:tests>>,copy,true>
            $<TARGET_RUNTIME_DLLS:tests> $<TARGET_FILE_DIR:tests>
        COMMAND_EXPAND_LISTS)
endif()
