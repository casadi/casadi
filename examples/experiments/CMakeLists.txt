find_package(yaml-cpp)
if (yaml-cpp_FOUND)
    add_executable(experiment experiment.cpp)
    target_link_libraries(experiment 
        PRIVATE 
            yaml-cpp
            cutest-loader
            panoc-alm::panoc-alm
    )

    add_executable(secondorderpanoc secondorderpanoc.cpp)
    target_link_libraries(secondorderpanoc 
        PRIVATE 
            yaml-cpp
            panoc-alm::casadi-loader
            panoc-alm::panoc-alm
    )
    casadi_function_codegen_python("himmelblau_functions" "codegen-himmelblau.py")
    add_dependencies(secondorderpanoc himmelblau_functions)

    add_executable(panocexperiment panocexperiment.cpp)
    target_link_libraries(panocexperiment 
        PRIVATE 
            yaml-cpp
            panoc-alm::panoc-alm
    )
endif()