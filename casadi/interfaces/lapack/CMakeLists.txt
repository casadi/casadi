cmake_minimum_required(VERSION 2.8.6)

# ifort on Windows
if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel" AND CMAKE_HOST_WIN32)
	add_definitions(-Ddgetrf_=DGETRF)
	add_definitions(-Ddgetrs_=DGETRS)
	add_definitions(-Ddgeequ_=DGEEQU)
	add_definitions(-Ddlaqge_=DLAQGE)
	add_definitions(-Ddgeqrf_=DGEQRF)
	add_definitions(-Ddormqr_=DORMQR)
	add_definitions(-Ddtrsm_=DTRSM)
endif()

# Linear solver using dense LU factorization in LAPACK
casadi_plugin(Linsol lapacklu
  lapack_lu.hpp
  lapack_lu.cpp
  lapack_lu_meta.cpp)
casadi_plugin_link_libraries(Linsol lapacklu ${LAPACK_LIBRARIES})

# Linear solver using dense QR factorization in LAPACK
casadi_plugin(Linsol lapackqr
  lapack_qr.hpp
  lapack_qr.cpp
  lapack_qr_meta.cpp
  )
casadi_plugin_link_libraries(Linsol lapackqr ${LAPACK_LIBRARIES})
