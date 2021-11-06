if (TARGET casadi)
    add_subdirectory(CasADi)
endif()

if (TARGET alpaqa::cutest-loader)
    add_subdirectory(CUTEst)
endif()

add_subdirectory(Presentation)
