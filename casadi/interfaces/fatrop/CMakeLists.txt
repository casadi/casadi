cmake_minimum_required(VERSION 3.10.2)

casadi_plugin(Conic fatrop
  fatrop_conic_interface.hpp
  fatrop_conic_interface.cpp
  fatrop_conic_interface_meta.cpp
  "${CMAKE_CURRENT_BINARY_DIR}/fatrop_conic_runtime_str.h"
  )

casadi_plugin_link_libraries(Conic fatrop fatrop::fatrop)

casadi_plugin(Nlpsol fatrop
  fatrop_interface.hpp
  fatrop_interface.cpp
  fatrop_interface_meta.cpp
  "${CMAKE_CURRENT_BINARY_DIR}/fatrop_runtime_str.h"
  )

casadi_plugin_link_libraries(Nlpsol fatrop fatrop::fatrop)


#add_executable(fatrop_test fatrop_test.cpp)
#target_link_libraries(fatrop_test fatrop casadi)


# Add the runtime sources to internal (installed separately)
set(RUNTIME_SRC "fatrop_conic_runtime.hpp")

string (REPLACE ";" "$<SEMICOLON>" ESCAPED_RUNTIME_SRC "${RUNTIME_SRC}")

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/fatrop_conic_runtime_str.h"
  COMMAND  ${CMAKE_COMMAND} -D OUTPUT="${CMAKE_CURRENT_BINARY_DIR}/fatrop_conic_runtime_str.h" -D SOURCES="${RUNTIME_SRC}" -P "${CMAKE_CURRENT_SOURCE_DIR}/../../generate_runtime.cmake"
  DEPENDS ${RUNTIME_SRC}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

# Add the runtime sources to internal (installed separately)
set(RUNTIME_SRC "fatrop_runtime.hpp")

string (REPLACE ";" "$<SEMICOLON>" ESCAPED_RUNTIME_SRC "${RUNTIME_SRC}")

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/fatrop_runtime_str.h"
  COMMAND  ${CMAKE_COMMAND} -D OUTPUT="${CMAKE_CURRENT_BINARY_DIR}/fatrop_runtime_str.h" -D SOURCES="${RUNTIME_SRC}" -P "${CMAKE_CURRENT_SOURCE_DIR}/../../generate_runtime.cmake"
  DEPENDS ${RUNTIME_SRC}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

target_include_directories(casadi_conic_fatrop
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)

target_include_directories(casadi_nlpsol_fatrop
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)