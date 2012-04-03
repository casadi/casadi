if(WITH_SUPERLU)
find_package (BLAS)
if(BLAS_FOUND)
add_subdirectory(superlu)
endif(BLAS_FOUND)
endif(WITH_SUPERLU)

if(LAPACK_FOUND)
add_subdirectory(lapack)
endif(LAPACK_FOUND)

if(SUNDIALS_FOUND)
add_subdirectory(sundials)
endif(SUNDIALS_FOUND)

IF (IPOPT_FOUND)
add_subdirectory(ipopt)
ENDIF (IPOPT_FOUND)

if(ACADO_FOUND)
add_subdirectory(acado)
add_subdirectory(qpoases)
endif(ACADO_FOUND)

IF (KNITRO_FOUND)
add_subdirectory(knitro)
ENDIF (KNITRO_FOUND)

IF (LIFTOPT_FOUND)
add_subdirectory(liftopt)
ENDIF (LIFTOPT_FOUND)

if (CPLEX_FOUND)
add_subdirectory(cplex)
endif (CPLEX_FOUND)

add_subdirectory(csparse)


if(WITH_GSL)
find_package (GSL)
find_package (BLAS)
if(GSL_FOUND)
if(BLAS_FOUND)
add_subdirectory(gsl)
endif(BLAS_FOUND)
endif(GSL_FOUND)
endif(WITH_GSL)

if(WITH_OOQP)
find_package (OOQP)
if(OOQP_FOUND)
add_subdirectory(ooqp)
endif(OOQP_FOUND)
endif(WITH_OOQP)

if(WITH_WORHP)
find_package (WORHP)
if(WORHP_FOUND)
add_subdirectory(worhp)
endif(WORHP_FOUND)
endif(WITH_WORHP)
