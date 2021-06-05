if (TARGET casadi)
    add_subdirectory(CasADi)
endif()

if (TARGET panoc-alm::cutest-loader)
    add_subdirectory(CUTEst)
endif()
