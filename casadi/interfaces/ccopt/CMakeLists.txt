cmake_minimum_required(VERSION 3.10.2)

casadi_plugin(Nlpsol madmpec
  madmpec_interface.hpp
  madmpec_interface.cpp
  madmpec_interface_meta.cpp
  "${CMAKE_CURRENT_BINARY_DIR}/madmpec_runtime_str.h"
  )

casadi_plugin_link_libraries(Nlpsol madmpec libmad::libMad)

#add_executable(madmpec_test madmpec_test.cpp)
#target_link_libraries(madmpec_test madmpec casadi)

# Add the runtime sources to internal (installed separately)
set(RUNTIME_SRC "madmpec_runtime.hpp")

string (REPLACE ";" "$<SEMICOLON>" ESCAPED_RUNTIME_SRC "${RUNTIME_SRC}")

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/madmpec_runtime_str.h"
  COMMAND  ${CMAKE_COMMAND} -D OUTPUT="${CMAKE_CURRENT_BINARY_DIR}/madmpec_runtime_str.h" -D SOURCES="${RUNTIME_SRC}" -P "${CMAKE_CURRENT_SOURCE_DIR}/../../generate_runtime.cmake"
  DEPENDS ${RUNTIME_SRC}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_directories(casadi_nlpsol_madmpec 
  PUBLIC
  ${LIBMAD_LIB_DIR}
)

#target_link_libraries(casadi_nlpsol_madmpec cmadmpec madmpec julia stdc++ dl m)

target_include_directories(casadi_nlpsol_madmpec
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)
include(CMakePrintHelpers)

message("CMAKE_CURRENT_BINARY_DIR = ${CMAKE_CURRENT_BINARY_DIR}")
