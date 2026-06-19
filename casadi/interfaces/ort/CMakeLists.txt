cmake_minimum_required(VERSION 3.10.2)

casadi_plugin(Onnx ort
  ort_interface.hpp
  ort_interface.cpp
  ort_runtime.h)

# The runtime header is included by generated C; install it for consumers
install(FILES ort_runtime.h
  DESTINATION ${INCLUDE_PREFIX}/casadi/interfaces/ort)

casadi_plugin_link_libraries(Onnx ort onnxruntime::onnxruntime)

target_include_directories(casadi_onnx_ort PRIVATE ${ONNX_RUNTIME_INCLUDE_DIR})
