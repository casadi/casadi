cmake_minimum_required(VERSION 3.10.2)

casadi_plugin(Nlpsol ccopt
  ccopt_interface.hpp
  ccopt_interface.cpp
  ccopt_interface_meta.cpp
  "${CMAKE_CURRENT_BINARY_DIR}/ccopt_runtime_str.h"
  )

casadi_plugin_link_libraries(Nlpsol ccopt libmad::libMad)

#add_executable(ccopt_test ccopt_test.cpp)
#target_link_libraries(ccopt_test ccopt casadi)

# Add the runtime sources to internal (installed separately)
set(RUNTIME_SRC "ccopt_runtime.hpp")

string (REPLACE ";" "$<SEMICOLON>" ESCAPED_RUNTIME_SRC "${RUNTIME_SRC}")

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/ccopt_runtime_str.h"
  COMMAND  ${CMAKE_COMMAND} -D OUTPUT="${CMAKE_CURRENT_BINARY_DIR}/ccopt_runtime_str.h" -D SOURCES="${RUNTIME_SRC}" -P "${CMAKE_CURRENT_SOURCE_DIR}/../../generate_runtime.cmake"
  DEPENDS ${RUNTIME_SRC}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_directories(casadi_nlpsol_ccopt 
  PUBLIC
  ${LIBMAD_LIB_DIR}
)

#target_link_libraries(casadi_nlpsol_ccopt cccopt ccopt julia stdc++ dl m)

target_include_directories(casadi_nlpsol_ccopt
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)
include(CMakePrintHelpers)

message("CMAKE_CURRENT_BINARY_DIR = ${CMAKE_CURRENT_BINARY_DIR}")
