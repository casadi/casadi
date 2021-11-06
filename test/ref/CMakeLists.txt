add_library(alpaqa-ref
    "src/panoc-ref.cpp"
    "src/fd.cpp"
)
target_include_directories(alpaqa-ref
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(alpaqa-ref PUBLIC alpaqa::alpaqa)