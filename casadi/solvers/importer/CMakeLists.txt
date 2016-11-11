cmake_minimum_required(VERSION 2.8.6)

if(WITH_DL AND NOT WIN32)
  # Simple just-in-time compiler, using shell commands
  casadi_plugin(Importer shell
    shell_compiler.hpp
    shell_compiler.cpp
    shell_compiler_meta.cpp)
endif()
