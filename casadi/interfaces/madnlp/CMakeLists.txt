cmake_minimum_required(VERSION 3.10.2)

casadi_plugin(Nlpsol madnlp
  madnlp_interface.hpp
  madnlp_interface.cpp
  madnlp_interface_meta.cpp
  "${CMAKE_CURRENT_BINARY_DIR}/madnlp_runtime_str.h"
  )

casadi_plugin_link_libraries(Nlpsol madnlp madnlp)


#add_executable(madnlp_test madnlp_test.cpp)
#target_link_libraries(madnlp_test madnlp casadi)

# Add the runtime sources to internal (installed separately)
set(RUNTIME_SRC "madnlp_runtime.hpp")

string (REPLACE ";" "$<SEMICOLON>" ESCAPED_RUNTIME_SRC "${RUNTIME_SRC}")

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/madnlp_runtime_str.h"
  COMMAND  ${CMAKE_COMMAND} -D OUTPUT="${CMAKE_CURRENT_BINARY_DIR}/madnlp_runtime_str.h" -D SOURCES="${RUNTIME_SRC}" -P "${CMAKE_CURRENT_SOURCE_DIR}/../../generate_runtime.cmake"
  DEPENDS ${RUNTIME_SRC}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

set(JULIA_INCLUDE_DIR /home/tom/.julia/juliaup/julia-1.10.4+0.x64.linux.gnu/include/julia/)
set(MADNLP_INCLUDE_DIR /home/tom/sources/_all/julia_casadi_interface/libmadnlp/cinterface/)
set(CMADNLP_INCLUDE_DIR /home/tom/sources/_all/julia_casadi_interface/libmadnlp/target/include/)

set(CMADNLP_LIB_DIR /home/tom/sources/_all/julia_casadi_interface/libmadnlp/target/lib/)

target_link_directories(casadi_nlpsol_madnlp 
  PUBLIC
  ${CMADNLP_LIB_DIR}
) 

target_link_libraries(casadi_nlpsol_madnlp cmadnlp madnlp julia stdc++ dl m)

target_include_directories(casadi_nlpsol_madnlp
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    ${MADNLP_INCLUDE_DIR}
    ${JULIA_INCLUDE_DIR}
)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")

