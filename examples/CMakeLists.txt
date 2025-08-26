cmake_minimum_required(VERSION 3.21...4.1)
project(alpaqa-examples C CXX)

if (PROJECT_IS_TOP_LEVEL)
    find_package(alpaqa REQUIRED CONFIG OPTIONAL_COMPONENTS CasADi Dl Extra)
endif()

add_subdirectory(util)

add_custom_target(examples ALL)
function(alpaqa_register_example TGT)
    target_link_libraries(${TGT} PRIVATE alpaqa::example-util)
    add_dependencies(examples ${TGT})
endfunction()

add_subdirectory(C++)
add_subdirectory(problems)