add_executable(casadi-rosenbrock main.cpp)
target_link_libraries(casadi-rosenbrock 
    PRIVATE
        alpaqa::casadi-loader
        alpaqa::alpaqa
        alpaqa::warnings)
casadi_function_codegen_python(rosenbrock_functions "codegen-rosenbrock.py")
target_compile_definitions(casadi-rosenbrock PRIVATE
    ROSENBROCK_FUNC_DLL=\"$<TARGET_FILE_NAME:rosenbrock_functions>\")
add_dependencies(casadi-rosenbrock rosenbrock_functions)