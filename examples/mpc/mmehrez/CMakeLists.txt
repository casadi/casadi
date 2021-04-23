if (TARGET casadi)
    add_executable(mmehrez mmehrez.cpp)
    target_link_libraries(mmehrez 
        PRIVATE
            panoc-alm::casadi-loader
            panoc-alm::panoc-alm)
    casadi_function_codegen_python("mmehrez_functions" "codegen-mmehrez.py")
    add_dependencies(mmehrez mmehrez_functions)
endif()