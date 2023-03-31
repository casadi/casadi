add_executable(custom-control-cpp-problem-example main.cpp)
target_link_libraries(custom-control-cpp-problem-example
    PRIVATE alpaqa::alpaqa alpaqa::warnings)