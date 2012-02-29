find_package(LLVM)
if(NOT LLVM_FOUND )
  message(FATAL_ERROR "LLVM package can't be found. Set CMAKE_PREFIX_PATH variable to LLVM's installation prefix.")
endif()

include_directories( ${LLVM_INCLUDE_DIR} )


# Simple example
add_executable(test_llvm test_llvm.cpp)
target_link_libraries(test_llvm ${LLVM_LIBRARIES})

