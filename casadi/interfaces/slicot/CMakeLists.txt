cmake_minimum_required(VERSION 2.6)

# the TRUE installs this to lib/casadi
casadi_library_or_plugin(casadi_slicot_interface TRUE
  slicot_tools.hpp
  slicot_tools.cpp)
target_link_libraries(casadi_slicot_interface ${SLICOT_LIBRARIES})
set_property(GLOBAL APPEND PROPERTY CASADI_SWIG_MODULES casadi_slicot_interface)

casadi_plugin(DpleSolver slicot
  psd_indef_dple_internal.hpp
  psd_indef_dple_internal.cpp
  psd_indef_dple_internal_meta.cpp)
casadi_plugin_link_libraries(DpleSolver slicot casadi_slicot_interface ${SLICOT_LIBRARIES})

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/slicot
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)

