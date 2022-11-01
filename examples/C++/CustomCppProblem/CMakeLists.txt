add_executable(custom-cpp-problem-example main.cpp)
target_link_libraries(custom-cpp-problem-example
    PRIVATE alpaqa::alpaqa alpaqa::warnings)