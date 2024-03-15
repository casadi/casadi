cmake_minimum_required(VERSION 3.10.2)


casadi_plugin(Conic daqp
  daqp_interface.hpp
  daqp_interface.cpp
  daqp_interface_meta.cpp
  "${CMAKE_CURRENT_BINARY_DIR}/daqp_runtime_str.h"
  )

casadi_plugin_link_libraries(Conic daqp daqp)

# Add the runtime sources to internal (installed separately)
set(RUNTIME_SRC "daqp_runtime.hpp")

string (REPLACE ";" "$<SEMICOLON>" ESCAPED_RUNTIME_SRC "${RUNTIME_SRC}")

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/daqp_runtime_str.h"
  COMMAND  ${CMAKE_COMMAND} -D OUTPUT="${CMAKE_CURRENT_BINARY_DIR}/daqp_runtime_str.h" -D SOURCES="${RUNTIME_SRC}" -P "${CMAKE_CURRENT_SOURCE_DIR}/../../generate_runtime.cmake"
  DEPENDS ${RUNTIME_SRC}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

target_include_directories(casadi_conic_daqp
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)