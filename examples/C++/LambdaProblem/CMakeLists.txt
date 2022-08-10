add_executable(lambda-problem-example main.cpp)
target_link_libraries(lambda-problem-example
    PRIVATE alpaqa::alpaqa alpaqa::warnings)