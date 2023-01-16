add_library(problem-c SHARED problem.c)
set_target_properties(problem-c PROPERTIES PREFIX "" DEBUG_POSTFIX "")

# DLL import/export
include(GenerateExportHeader)
generate_export_header(problem-c
    EXPORT_FILE_NAME problem-c-export.h)
target_include_directories(problem-c PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
target_link_libraries(problem-c PRIVATE alpaqa::dl-headers alpaqa::warnings)

add_executable(dl-problem-example main.cpp)
target_link_libraries(dl-problem-example
    PRIVATE alpaqa::alpaqa alpaqa::warnings alpaqa::dl-loader)
target_compile_definitions(dl-problem-example PRIVATE
    DLPROBLEM_DLL=\"$<TARGET_FILE_NAME:problem-c>\")
add_dependencies(dl-problem-example problem-c)

add_executable(test-matmul "test-matmul.c")
target_link_libraries(test-matmul PRIVATE alpaqa::dl-headers alpaqa::warnings)
