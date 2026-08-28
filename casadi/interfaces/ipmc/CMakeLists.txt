cmake_minimum_required(VERSION 3.16.3)

casadi_plugin(Nlpsol ipmc
  ipmc_interface.hpp
  ipmc_interface.cpp
  ipmc_interface_meta.cpp
  "${CMAKE_CURRENT_BINARY_DIR}/ipmc_runtime_str.h"
  )

casadi_plugin_link_libraries(Nlpsol ipmc ipmc::ipmc blasfeo)

# Add the runtime sources to internal (installed separately)
set(RUNTIME_SRC "ipmc_runtime.hpp")

string (REPLACE ";" "$<SEMICOLON>" ESCAPED_RUNTIME_SRC "${RUNTIME_SRC}")

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/ipmc_runtime_str.h"
  COMMAND  ${CMAKE_COMMAND} -D OUTPUT="${CMAKE_CURRENT_BINARY_DIR}/ipmc_runtime_str.h" -D SOURCES="${RUNTIME_SRC}" -P "${CMAKE_CURRENT_SOURCE_DIR}/../../generate_runtime.cmake"
  DEPENDS ${RUNTIME_SRC}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

target_include_directories(casadi_nlpsol_ipmc
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)

if(WITH_MOCKUP_IPMC)
  set(_ipmc_libdir "")
else()
  get_filename_component(_ipmc_libdir "${IPMC_LIBRARY}" DIRECTORY)
endif()

set_target_properties(casadi_nlpsol_ipmc PROPERTIES
  BUILD_RPATH "${_ipmc_libdir}"
  INSTALL_RPATH "${_ipmc_libdir};${CMAKE_INSTALL_RPATH}"
  INSTALL_RPATH_USE_LINK_PATH FALSE
  BUILD_WITH_INSTALL_RPATH FALSE
)
