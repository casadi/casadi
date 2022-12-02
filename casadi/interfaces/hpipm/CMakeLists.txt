cmake_minimum_required(VERSION 2.8.6)

casadi_plugin(Conic hpipm
  hpipm_interface.hpp
  hpipm_interface.cpp
  hpipm_interface_meta.cpp
  "${CMAKE_CURRENT_BINARY_DIR}/hpipm_runtime_str.h"
  )

casadi_plugin_link_libraries(Conic hpipm hpipm)

# Add the runtime sources to internal (installed separately)
set(RUNTIME_SRC "hpipm_runtime.hpp")

string (REPLACE ";" "$<SEMICOLON>" ESCAPED_RUNTIME_SRC "${RUNTIME_SRC}")

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/hpipm_runtime_str.h"
  COMMAND  ${CMAKE_COMMAND} -D OUTPUT="${CMAKE_CURRENT_BINARY_DIR}/hpipm_runtime_str.h" -D SOURCES="${RUNTIME_SRC}" -P "${CMAKE_CURRENT_SOURCE_DIR}/../../generate_runtime.cmake"
  DEPENDS ${RUNTIME_SRC}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

target_include_directories(casadi_conic_hpipm
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)