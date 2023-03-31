add_executable(custom-control-cpp-problem main.cpp)
target_link_libraries(custom-control-cpp-problem
    PRIVATE alpaqa::alpaqa alpaqa::warnings)