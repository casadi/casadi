cmake_minimum_required(VERSION 2.8.6)

# MATLAB front-end
add_custom_target(matlab DEPENDS casadi_matlab)

# Find packages
find_package(MATLAB REQUIRED)

include_directories(${MATLAB_INCLUDE_DIR})
add_definitions(${MATLAB_EXTRA_CXXFLAGS})

# Generate SWIG wrapper
set_source_files_properties(../casadi.i  PROPERTIES  CPLUSPLUS ON)
swig_add_module(casadi_matlab matlab ../casadi.i)
