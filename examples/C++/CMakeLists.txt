add_subdirectory(LambdaProblem)

include(CheckLanguage)
check_language(Fortran)
if (CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
    add_subdirectory(FortranProblem)
endif()

find_package(casadi)
if (TARGET casadi)
    add_subdirectory(CasADi)
endif()