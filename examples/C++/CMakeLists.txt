add_subdirectory(CustomCppProblem)

include(CheckLanguage)
check_language(Fortran)
if (CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
    add_subdirectory(FortranProblem)
endif()

find_package(casadi)
if (TARGET casadi AND TARGET alpaqa::casadi-loader)
    add_subdirectory(CasADi)
endif()
