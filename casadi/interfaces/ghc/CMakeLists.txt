cmake_minimum_required(VERSION 3.16.3)

casadi_plugin(Filesystem ghc
  ghc.hpp
  ghc.cpp
  ghc_meta.cpp)
casadi_plugin_link_libraries(Filesystem ghc ghcFilesystem::ghc_filesystem)