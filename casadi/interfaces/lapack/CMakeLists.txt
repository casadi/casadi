cmake_minimum_required(VERSION 2.8.6)
include_directories(../)

# Linear solver using dense LU factorization in LAPACK
casadi_plugin(Linsol lapacklu
  lapack_lu_dense.hpp
  lapack_lu_dense.cpp
  lapack_lu_dense_meta.cpp)
casadi_plugin_link_libraries(Linsol lapacklu ${LAPACK_LIBRARIES})

# Linear solver using dense QR factorization in LAPACK
casadi_plugin(Linsol lapackqr
  lapack_qr_dense.hpp
  lapack_qr_dense.cpp
  lapack_qr_dense_meta.cpp
  )
casadi_plugin_link_libraries(Linsol lapackqr ${LAPACK_LIBRARIES})
