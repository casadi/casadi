add_executable(casadi-rosenbrock main.cpp)
target_link_libraries(casadi-rosenbrock PRIVATE casadi panoc-alm::panoc-alm)

find_package(Python3 REQUIRED COMPONENTS Interpreter)

set(PY_CODE_GEN_FILE "${CMAKE_CURRENT_SOURCE_DIR}/codegen-rosenbrock.py")
add_custom_command(OUTPUT rosenbrock_functions.c
                   COMMAND Python3::Interpreter "${PY_CODE_GEN_FILE}"
                   MAIN_DEPENDENCY "${PY_CODE_GEN_FILE}")
add_library(rosenbrock_functions SHARED rosenbrock_functions.c)
set_target_properties(rosenbrock_functions PROPERTIES DEBUG_POSTFIX ""
                                                      ASAN_POSTFIX ""
                                                      TSAN_POSTFIX "")
target_compile_options(rosenbrock_functions PRIVATE "-Wno-unused-parameter")
add_dependencies(casadi-rosenbrock rosenbrock_functions)