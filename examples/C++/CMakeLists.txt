add_subdirectory(CustomCppProblem)
add_subdirectory(SimpleUnconstrProblem)
add_subdirectory(Advanced)

include(CheckLanguage)
check_language(Fortran)
if (CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
    add_subdirectory(FortranProblem)
endif()

find_package(casadi)
if (TARGET casadi AND NOT TARGET casadi::casadi)
    add_library(casadi::casadi ALIAS casadi)  # Support CasADi <3.7
endif()
if (TARGET casadi::casadi AND TARGET alpaqa::casadi-loader)
    add_subdirectory(CasADi)
endif()

if (TARGET alpaqa::dl-loader)
    add_subdirectory(DLProblem)
endif()

if (ALPAQA_WITH_OCP)
    add_subdirectory(CustomControlCppProblem)
endif()

option(ALPAQA_WITH_CUTEST_EXAMPLES
    "Build the CUTEst example problems" Off)
if (ALPAQA_WITH_CUTEST_EXAMPLES)
    add_subdirectory(CUTEst)
endif()
