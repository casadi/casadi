cmake_minimum_required(VERSION 2.6)

include_directories(../)
include_directories(${SNOPT_INCLUDE_DIR})

casadi_plugin(NlpSolver snopt
  snopt_interface.hpp
  snopt_interface.cpp
  wsnopt.hpp
  snopt_interface_meta.cpp)

if(WITH_SNOPT_FORCE_DUMMY)
  if(ENABLE_STATIC)
  add_library(snopt7 STATIC empty_lib.cpp)
  endif()
  if(ENABLE_SHARED)
  add_library(snopt7 SHARED empty_lib.cpp)
  endif()
  casadi_plugin_link_libraries(NlpSolver snopt snopt7)
else()
  casadi_plugin_link_libraries(NlpSolver snopt ${SNOPT_LIBRARIES})
endif()

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/snopt
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
