if(LAPACK_FOUND)
  add_subdirectory(lapack)
endif()

if(WITH_SUNDIALS)
  add_subdirectory(sundials)
endif()

if(IPOPT_FOUND)
  add_subdirectory(ipopt)
endif()

if(ACADO_FOUND)
  add_subdirectory(acado)
endif(ACADO_FOUND)

if(QPOASES_FOUND)
  add_subdirectory(qpoases)
endif()

if(KNITRO_FOUND)
  add_subdirectory(knitro)
endif()

if(LIFTOPT_FOUND)
  add_subdirectory(liftopt)
endif()

if(CPLEX_FOUND)
  add_subdirectory(cplex)
endif()

if(WITH_CSPARSE)
  add_subdirectory(csparse)
endif()

if(WITH_DSDP)
  add_subdirectory(dsdp)
endif()

# if(GSL_FOUND)
#   add_subdirectory(gsl)
# endif()

if(OOQP_FOUND)
  add_subdirectory(ooqp)
endif()

if(WORHP_FOUND)
  add_subdirectory(worhp)
endif()
