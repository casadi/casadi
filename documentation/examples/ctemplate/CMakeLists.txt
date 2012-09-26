if(ENABLE_SHARED)
  set(EXTRA_LIBS ${CMAKE_DL_LIBS})
else()
  set(EXTRA_LIBS "")
endif()

if(LAPACK_FOUND AND BLAS_FOUND)
add_executable(ctemplate ctemplate.cpp)
target_link_libraries(ctemplate
  sundials_interface 
  lapack_interface
  csparse_interface
  casadi
  ${SUNDIALS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  csparse
  ${BLAS_LIBRARIES}
  ${EXTRA_LIBS}
)
endif()
