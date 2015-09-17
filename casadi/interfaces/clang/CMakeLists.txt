cmake_minimum_required(VERSION 2.8.6)
include_directories(${CLANG_INCLUDE_DIR})
add_definitions(${LLVM_DEFINITIONS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CLANG_CXX_FLAGS}")

casadi_plugin(JitCompiler clang
  clang_interface.hpp
  clang_interface.cpp
  clang_interface_meta.cpp)

casadi_plugin_link_libraries(JitCompiler clang ${CLANG_LIBRARIES})

# Install resource directory
INSTALL( DIRECTORY ${LLVM_INSTALL_PREFIX}/lib/clang DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/casadi/jit/)

# Install libc++
INSTALL( DIRECTORY ${LLVM_INSTALL_PREFIX}/include/c++ DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/casadi/jit/)

