find_package(IPOPT QUIET)
find_package(KNITRO QUIET)
find_package(BLAS QUIET)
find_package(LAPACK QUIET)
find_package(LIFTOPT QUIET)

if(ENABLE_SHARED)
  set(EXTRA_LIBS ${CMAKE_DL_LIBS})
else()
  set(EXTRA_LIBS "")
endif()

if(WITH_SUNDIALS AND LAPACK_FOUND AND WITH_SUPERLU AND BLAS_FOUND)
add_executable(ctemplate ctemplate.cpp)
target_link_libraries(ctemplate
  sundials_interface
  superlu_interface 
  lapack_interface
  csparse_interface
  casadi
  ${SUNDIALS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  superlu
  csparse
  ${BLAS_LIBRARIES}
  ${EXTRA_LIBS}
)
endif()
