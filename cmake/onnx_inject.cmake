# Injected into ONNX build via CMAKE_PROJECT_INCLUDE
# Forces ONNX to use our pre-configured protobuf

find_package(Protobuf CONFIG REQUIRED)

message(STATUS "ONNX inject: Found protobuf via CONFIG")
message(STATUS "ONNX inject: protobuf::libprotobuf = $<TARGET_FILE:protobuf::libprotobuf>")

# If a custom protoc executable was specified, make it available
if(DEFINED ONNX_CUSTOM_PROTOC_EXECUTABLE)
  message(STATUS "ONNX inject: Using custom protoc: ${ONNX_CUSTOM_PROTOC_EXECUTABLE}")
endif()
