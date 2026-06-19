cmake_minimum_required(VERSION 3.16.3)

casadi_plugin(Nlpsol uno
  uno_interface.hpp
  uno_interface.cpp
  uno_interface_meta.cpp
  uno_runtime.hpp
  "${CMAKE_CURRENT_BINARY_DIR}/uno_runtime_str.h")

casadi_plugin_link_libraries(Nlpsol uno uno)

# Stringify uno_runtime.hpp -> uno_runtime_str.h, mirroring fatrop's pattern.
set(RUNTIME_SRC "uno_runtime.hpp")
string(REPLACE ";" "$<SEMICOLON>" ESCAPED_RUNTIME_SRC "${RUNTIME_SRC}")
add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/uno_runtime_str.h"
  COMMAND ${CMAKE_COMMAND}
    -D OUTPUT="${CMAKE_CURRENT_BINARY_DIR}/uno_runtime_str.h"
    -D SOURCES="${RUNTIME_SRC}"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/../../generate_runtime.cmake"
  DEPENDS ${RUNTIME_SRC}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

target_include_directories(casadi_nlpsol_uno
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
