cmake_minimum_required(VERSION 2.6)

set(SLICOT_INTERFACE_SRCS
  slicot_tools.hpp
  slicot_tools.cpp
  psd_indef_dple_solver.hpp
  psd_indef_dple_solver.cpp
  psd_indef_dple_internal.hpp
  psd_indef_dple_internal.cpp
)

if(ENABLE_STATIC)
add_library(casadi_slicot_interface STATIC ${SLICOT_INTERFACE_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(casadi_slicot_interface SHARED ${SLICOT_INTERFACE_SRCS})
endif(ENABLE_SHARED)
install(TARGETS casadi_slicot_interface
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY ./ 
   DESTINATION include/casadi/interfaces/slicot
   FILES_MATCHING PATTERN "*.hpp"
   PATTERN ".svn" EXCLUDE
)
