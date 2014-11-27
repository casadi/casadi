cmake_minimum_required(VERSION 2.8.6)

casadi_plugin(DpleSolver slicot
  psd_indef_dple_internal.hpp
  psd_indef_dple_internal.cpp
  psd_indef_dple_internal_meta.cpp)
casadi_plugin_link_libraries(DpleSolver slicot ${SLICOT_LIBRARIES})
