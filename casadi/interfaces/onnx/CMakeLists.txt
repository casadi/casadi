cmake_minimum_required(VERSION 3.10.2)

# ONNX Translator plugin
if(WITH_ONNX)
  casadi_plugin(Translator onnx
    onnx_translator.hpp
    onnx_translator.cpp
    onnx_translator_meta.cpp
    onnx_utilities.cpp
    onnx_operations.cpp
    onnx_import.cpp
    onnx_export.cpp
  )

  # Link against ONNX libraries
  casadi_plugin_link_libraries(Translator onnx ONNX::onnx ONNX::onnx_proto)

  # Force C++17 for ONNX compatibility
  set_target_properties(casadi_translator_onnx PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
  )
endif()
