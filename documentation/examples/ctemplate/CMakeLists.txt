if(ENABLE_SHARED)
  set(EXTRA_LIBS ${CMAKE_DL_LIBS})
else()
  set(EXTRA_LIBS "")
endif()

if(LAPACK_FOUND AND BLAS_FOUND AND WITH_CSPARSE)
add_executable(ctemplate ctemplate.cpp)
target_link_libraries(ctemplate
  casadi_sundials_interface 
  casadi_lapack_interface
  casadi_csparse_interface
  casadi
  ${SUNDIALS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${CSPARSE_LIBRARIES}
  ${BLAS_LIBRARIES}
  ${EXTRA_LIBS}
  ${CASADI_DEPENDENCIES}
)
endif()
