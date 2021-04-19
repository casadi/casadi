find_package(yaml-cpp)
if (yaml-cpp_FOUND)
    add_executable(poster poster.cpp)
    target_link_libraries(poster
        PRIVATE 
            yaml-cpp
            panoc-alm::panoc-alm
    )
    add_executable(pga-vs-panoc-vs-newton pga-vs-panoc-vs-newton.cpp)
    target_link_libraries(pga-vs-panoc-vs-newton
        PRIVATE 
            yaml-cpp
            panoc-alm::panoc-alm
    )
endif()