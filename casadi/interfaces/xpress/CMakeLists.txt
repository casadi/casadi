cmake_minimum_required(VERSION 3.10.2)

include_directories(../)

casadi_plugin(Conic xpress
  xpress_interface.hpp
  xpress_interface.cpp
  xpress_interface_meta.cpp
  xpress_runtime.hpp
  "${CMAKE_CURRENT_BINARY_DIR}/xpress_runtime_str.h"
  )

casadi_plugin_link_libraries(Conic xpress xpress::xpress)

# Generate runtime string header for codegen
set(RUNTIME_SRC "xpress_runtime.hpp")

string (REPLACE ";" "$<SEMICOLON>" ESCAPED_RUNTIME_SRC "${RUNTIME_SRC}")

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/xpress_runtime_str.h"
  COMMAND  ${CMAKE_COMMAND} -D OUTPUT="${CMAKE_CURRENT_BINARY_DIR}/xpress_runtime_str.h" -D SOURCES="${RUNTIME_SRC}" -P "${CMAKE_CURRENT_SOURCE_DIR}/../../generate_runtime.cmake"
  DEPENDS ${RUNTIME_SRC}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

target_include_directories(casadi_conic_xpress
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)
