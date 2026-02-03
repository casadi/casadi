cmake_minimum_required(VERSION 3.16.3)

casadi_plugin(Archiver libzip
  libzip.hpp
  libzip.cpp
  libzip_meta.cpp)
casadi_plugin_link_libraries(Archiver libzip libzip::zip ZLIB::libz)
