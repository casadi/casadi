cmake_minimum_required(VERSION 3.10.2)

# ONNX GraphModel backend: the protobuf graph engine, exposed via the GraphModel family
if(WITH_ONNX)
  casadi_plugin(GraphModel onnx
    onnx_model.hpp
    onnx_model.cpp
    onnx_utilities.cpp
    onnx_operations.cpp
    onnx_import.cpp
    onnx_export.cpp
  )

  # Link against ONNX libraries
  casadi_plugin_link_libraries(GraphModel onnx ONNX::onnx ONNX::onnx_proto)

  # Force C++17 for ONNX compatibility
  set_target_properties(casadi_graphmodel_onnx PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
  )
endif()
