find_package(yaml-cpp)
if (yaml-cpp_FOUND)
    add_executable(experiment experiment.cpp)
    target_link_libraries(experiment 
        PRIVATE 
            yaml-cpp
            cutest-loader
            panoc-alm::panoc-alm
    )

    add_executable(panocexperiment panocexperiment.cpp)
    target_link_libraries(panocexperiment 
        PRIVATE 
            yaml-cpp
            panoc-alm::panoc-alm
    )
endif()