cmake_minimum_required(VERSION 2.8.6)

# MATLAB front-end
add_custom_target(matlab DEPENDS casadiMATLAB_wrap)

# Find MATLAB
find_package(MATLAB REQUIRED)

include_directories(${MATLAB_INCLUDE_DIR})
add_definitions(${MATLAB_EXTRA_CXXFLAGS})

# Generate SWIG wrapper
set_source_files_properties(../casadi.i  PROPERTIES  CPLUSPLUS ON)
swig_add_module(casadiMATLAB_wrap matlab ../casadi.i)
swig_link_libraries(casadiMATLAB_wrap ${MATLAB_LIBRARIES} casadi)
set_target_properties(casadiMATLAB_wrap PROPERTIES PREFIX "" SUFFIX .${MATLAB_MEX_EXT})

# Install the generated front-end to ${CMAKE_INSTALL_PREFIX}/matlab_tb
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/+casadi DESTINATION ${CMAKE_INSTALL_PREFIX}/matlab_tb)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/SwigRef.m DESTINATION ${CMAKE_INSTALL_PREFIX}/matlab_tb)
install(TARGETS casadiMATLAB_wrap DESTINATION ${CMAKE_INSTALL_PREFIX}/matlab_tb)
