cmake_minimum_required(VERSION 2.6)

include_directories(../)
include_directories(${SNOPT_INCLUDE_DIR})

casadi_plugin(NlpSolver snopt
  snopt_interface.hpp
  snopt_interface.cpp
  snopt.h
  snoptProblem.hpp
  snopt_interface_meta.cpp)

if(WITH_SNOPT_FORCE_DUMMY)
  if(ENABLE_STATIC)
  add_library(snopt7 STATIC empty_lib.cpp)
  add_library(snopt7_cpp STATIC empty_lib.cpp)
  endif()
  if(ENABLE_SHARED)
  add_library(snopt7 SHARED empty_lib.cpp)
  add_library(snopt7_cpp SHARED empty_lib.cpp)
  endif()
  casadi_plugin_link_libraries(NlpSolver snopt snopt7 snopt7_cpp)
else()
  casadi_plugin_link_libraries(NlpSolver snopt ${SNOPT_LIBRARIES})
endif()

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/snopt
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
