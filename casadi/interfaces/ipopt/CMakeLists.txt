cmake_minimum_required(VERSION 3.10.2)

set(NLPSOL_IPOPT_SRCS
  ipopt_interface.hpp
  ipopt_interface.cpp
  ipopt_nlp.hpp
  ipopt_nlp.cpp
  ipopt_interface_meta.cpp
  "${CMAKE_CURRENT_BINARY_DIR}/ipopt_runtime_str.h")

if(WITH_IPOPT_CALLBACK)
  add_definitions(-DWITH_IPOPT_CALLBACK)
endif()

casadi_plugin(Nlpsol ipopt ${NLPSOL_IPOPT_SRCS})
casadi_plugin_link_libraries(Nlpsol ipopt ipopt)

# Add the runtime sources to internal (installed separately)
set(RUNTIME_SRC "ipopt_runtime.hpp")

string (REPLACE ";" "$<SEMICOLON>" ESCAPED_RUNTIME_SRC "${RUNTIME_SRC}")

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/ipopt_runtime_str.h"
  COMMAND  ${CMAKE_COMMAND} -D OUTPUT="${CMAKE_CURRENT_BINARY_DIR}/ipopt_runtime_str.h" -D SOURCES="${RUNTIME_SRC}" -P "${CMAKE_CURRENT_SOURCE_DIR}/../../generate_runtime.cmake"
  DEPENDS ${RUNTIME_SRC}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

target_include_directories(casadi_nlpsol_ipopt
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)

#set_property(TARGET casadi_nlpsol_ipopt PROPERTY PKG_CONFIG_REQUIRES ipopt)
