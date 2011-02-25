if(WITH_SUPERLU)
find_package (BLAS)
if(BLAS_FOUND)
add_subdirectory(superlu)
endif(BLAS_FOUND)
endif(WITH_SUPERLU)

find_package (LAPACK)
if(LAPACK_FOUND)
add_subdirectory(lapack)
endif(LAPACK_FOUND)

find_package(SUNDIALS)
if(SUNDIALS_FOUND)
add_subdirectory(sundials)
endif(SUNDIALS_FOUND)

find_package(IPOPT)
IF (IPOPT_FOUND)
add_subdirectory(ipopt)
ENDIF (IPOPT_FOUND)

find_package(ACADO)
if(ACADO_FOUND)
add_subdirectory(acado)
endif(ACADO_FOUND)

find_package(KNITRO)
IF (KNITRO_FOUND)
add_subdirectory(knitro)
ENDIF (KNITRO_FOUND)

find_package(LIFTOPT)
IF (LIFTOPT_FOUND)
add_subdirectory(liftopt)
ENDIF (LIFTOPT_FOUND)

find_package(CPLEX)
if (CPLEX_FOUND)
add_subdirectory(cplex)
endif (CPLEX_FOUND)

add_subdirectory(csparse)
