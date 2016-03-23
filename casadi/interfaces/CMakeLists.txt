if(LAPACK_FOUND)
  add_subdirectory(lapack)
endif()

if(SUNDIALS_FOUND)
  add_subdirectory(sundials)
endif()

if(IPOPT_FOUND)
  add_subdirectory(ipopt)
endif()

if(QPOASES_FOUND)
  add_subdirectory(qpoases)
endif()

if(KNITRO_FOUND)
  add_subdirectory(knitro)
endif()

if(CPLEX_FOUND)
  add_subdirectory(cplex)
endif()

if(CSPARSE_FOUND)
  add_subdirectory(csparse)
endif()

if(CLANG_FOUND)
  add_subdirectory(clang)
endif()

if(OOQP_FOUND)
  add_subdirectory(ooqp)
endif()

if(GUROBI_FOUND)
  add_subdirectory(gurobi)
endif()

if(WORHP_FOUND)
  add_subdirectory(worhp)
endif()

if(SQIC_FOUND)
  add_subdirectory(sqic)
endif()

if(SNOPT_FOUND OR WITH_SNOPT_FORCE_DUMMY)
  add_subdirectory(snopt)
endif()

if(TINYXML_FOUND)
  add_subdirectory(tinyxml)
endif()

set(LINT_TARGETS ${LINT_TARGETS} PARENT_SCOPE)
set(SPELL_TARGETS ${SPELL_TARGETS} PARENT_SCOPE)
set(CASADI_MODULES ${CASADI_MODULES} PARENT_SCOPE)
