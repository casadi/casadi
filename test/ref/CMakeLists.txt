add_library(panoc-alm-ref
    "src/panoc-ref.cpp"
    "src/fd.cpp"
)
target_include_directories(panoc-alm-ref
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(panoc-alm-ref PUBLIC panoc-alm::panoc-alm)