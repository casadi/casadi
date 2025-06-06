if(WITH_LAPACK)
  add_subdirectory(lapack)
endif()

if(WITH_SUNDIALS)
  add_subdirectory(sundials)
endif()

if(WITH_SLEQP)
  add_subdirectory(sleqp)
endif()

if(WITH_IPOPT)
  add_subdirectory(ipopt)
endif()

if(WITH_MADNLP)
  add_subdirectory(madnlp)
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

if(WITH_CBC)
  add_subdirectory(cbc)
endif()

if(WITH_CSPARSE)
  add_subdirectory(csparse)
endif()

if(WITH_CLANG)
  add_subdirectory(clang)
endif()

if(WITH_HIGHS)
  add_subdirectory(highs)
endif()

if(WITH_CLARABEL)
  add_subdirectory(clarabel)
endif()

if(WITH_DAQP)
  add_subdirectory(daqp)
endif()

if(WITH_HSL)
  add_subdirectory(hsl)
endif()

if(WITH_MUMPS)
  add_subdirectory(mumps)
endif()

if(WITH_OOQP)
  add_subdirectory(ooqp)
endif()

if(WITH_GUROBI)
  add_subdirectory(gurobi)
endif()

if(WITH_WORHP)
  add_subdirectory(worhp)
endif()

if(WITH_SQIC)
  add_subdirectory(sqic)
endif()

if(WITH_SNOPT)
  add_subdirectory(snopt)
endif()

if(WITH_TINYXML)
  add_subdirectory(tinyxml)
endif()

if(WITH_BLOCKSQP)
  add_subdirectory(blocksqp)
endif()

if(WITH_HPMPC)
  add_subdirectory(hpmpc)
endif()

if(WITH_HPIPM)
  add_subdirectory(hpipm)
endif()

if(WITH_FATROP)
  add_subdirectory(fatrop)
endif()

if(WITH_SUPERSCS)
  add_subdirectory(superscs)
endif()

if(WITH_SLICOT)
  add_subdirectory(slicot)
endif()

if(WITH_AMPL)
  add_subdirectory(ampl)
endif()

if(WITH_PROXQP)
  add_subdirectory(proxqp)
endif()

if(WITH_OSQP)
  add_subdirectory(osqp)
endif()

if(WITH_ALPAQA)
  add_subdirectory(alpaqa)
endif()

if(WITH_MATLAB_IPC)
  add_subdirectory(matlab_ipc)
endif()

if(WITH_LIBZIP)
  add_subdirectory(libzip)
endif()

if(WITH_GHC_FILESYSTEM)
  add_subdirectory(ghc)
endif()

if(WITH_LINT)
  set(LINT_TARGETS ${LINT_TARGETS} PARENT_SCOPE)
endif()

if(WITH_SPELL)
  set(SPELL_TARGETS ${SPELL_TARGETS} PARENT_SCOPE)
endif()

set(CASADI_MODULES ${CASADI_MODULES} PARENT_SCOPE)
