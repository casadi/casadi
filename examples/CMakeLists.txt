if (TARGET casadi)
    add_subdirectory(CasADi)
endif()

if (TARGET panoc-alm::cutest-loader)
    add_subdirectory(CUTEst)
endif()

add_subdirectory(experiments)
add_subdirectory(poster)