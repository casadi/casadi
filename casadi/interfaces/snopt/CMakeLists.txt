cmake_minimum_required(VERSION 2.6)

include_directories(../)
include_directories(${SNOPT_INCLUDE_DIR})

set(SNOPT_INTERFACE_SRCS
  snopt_internal.hpp
  snopt_internal.cpp
  snopt_solver.hpp
  snopt_solver.cpp
  wsnopt.hpp)

casadi_library(casadi_snopt_interface ${SNOPT_INTERFACE_SRCS})

if(WITH_SNOPT_FORCE_DUMMY)
  if(ENABLE_STATIC)
  add_library(snopt7 STATIC empty_lib.cpp)
  endif()
  if(ENABLE_SHARED)
  add_library(snopt7 SHARED empty_lib.cpp)
  endif()
  target_link_libraries(casadi_snopt_interface ${CMAKE_DL_LIBS} snopt7)
else()
  target_link_libraries(casadi_snopt_interface ${CMAKE_DL_LIBS} ${SNOPT_LIBRARIES})
endif()

install(DIRECTORY ./
  DESTINATION include/casadi/interfaces/snopt
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
