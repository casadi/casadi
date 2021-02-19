add_executable(casadi-rosenbrock main.cpp)
target_link_libraries(casadi-rosenbrock PRIVATE casadi panoc-alm::panoc-alm)
casadi_function_codegen_python("rosenbrock_functions" "codegen-rosenbrock.py")
add_dependencies(casadi-rosenbrock rosenbrock_functions)