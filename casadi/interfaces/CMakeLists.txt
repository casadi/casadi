if(LAPACK_FOUND)
  add_subdirectory(lapack)
endif()

if(SUNDIALS_FOUND)
  add_subdirectory(sundials)
endif()

if(WITH_IPOPT)
  add_subdirectory(ipopt)
endif()

if(WITH_BONMIN)
  add_subdirectory(bonmin)
endif()

if(WITH_QPOASES)
  add_subdirectory(qpoases)
endif()

if(WITH_KNITRO)
  add_subdirectory(knitro)
endif()

if(WITH_CPLEX)
  add_subdirectory(cplex)
endif()

if(WITH_CLP)
  add_subdirectory(clp)
endif()

if(WITH_CSPARSE)
  add_subdirectory(csparse)
endif()

if(CLANG_FOUND)
  add_subdirectory(clang)
endif()

if(HSL_FOUND)
  add_subdirectory(hsl)
endif()

if(OOQP_FOUND)
  add_subdirectory(ooqp)
endif()

if(WITH_GUROBI)
  add_subdirectory(gurobi)
endif()

if(WITH_WORHP)
  add_subdirectory(worhp)
endif()

if(SQIC_FOUND)
  add_subdirectory(sqic)
endif()

if(WITH_SNOPT)
  add_subdirectory(snopt)
endif()

if(TINYXML_FOUND)
  add_subdirectory(tinyxml)
endif()

if(BLOCKSQP_FOUND)
  add_subdirectory(blocksqp)
endif()

if(WITH_HPMPC)
  add_subdirectory(hpmpc)
endif()

if(SLICOT_FOUND)
  add_subdirectory(slicot)
endif()

if(WITH_LINT)
  set(LINT_TARGETS ${LINT_TARGETS} PARENT_SCOPE)
endif()

if(WITH_SPELL)
  set(SPELL_TARGETS ${SPELL_TARGETS} PARENT_SCOPE)
endif()

set(CASADI_MODULES ${CASADI_MODULES} PARENT_SCOPE)
