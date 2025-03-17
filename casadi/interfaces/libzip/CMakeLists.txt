cmake_minimum_required(VERSION 3.10.2)

casadi_plugin(Archiver libzip
  libzip.hpp
  libzip.cpp
  libzip_meta.cpp)
casadi_plugin_link_libraries(Archiver libzip libzip::zip ZLIB::libz)
