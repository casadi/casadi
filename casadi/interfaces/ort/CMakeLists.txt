cmake_minimum_required(VERSION 3.10.2)

casadi_plugin(Onnx ort
  ort_interface.hpp
  ort_interface.cpp
  "${CMAKE_CURRENT_BINARY_DIR}/ort_runtime_str.h")

# Install the runtime source alongside the interface headers
install(FILES ort_runtime.hpp
  DESTINATION ${INCLUDE_PREFIX}/casadi/interfaces/ort)

casadi_plugin_link_libraries(Onnx ort onnxruntime::onnxruntime)

target_include_directories(casadi_onnx_ort PRIVATE ${ONNX_RUNTIME_INCLUDE_DIR})

# Add the runtime sources to internal (installed separately)
set(RUNTIME_SRC "ort_runtime.hpp")

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/ort_runtime_str.h"
  COMMAND ${CMAKE_COMMAND}
    -D OUTPUT="${CMAKE_CURRENT_BINARY_DIR}/ort_runtime_str.h"
    -D SOURCES="${RUNTIME_SRC}"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/../../generate_runtime.cmake"
  DEPENDS ${RUNTIME_SRC}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

target_include_directories(casadi_onnx_ort
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)
