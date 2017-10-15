cmake_minimum_required(VERSION 2.8.6)

casadi_plugin(Linsol symbolicqr
  symbolic_qr.hpp symbolic_qr.cpp symbolic_qr_meta.cpp
)

casadi_plugin(Linsol qr
  qr.hpp qr.cpp qr_meta.cpp
)

casadi_plugin(Linsol lsqr
  lsqr.hpp lsqr.cpp lsqr_meta.cpp
)
