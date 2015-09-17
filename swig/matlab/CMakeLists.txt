cmake_minimum_required(VERSION 2.8.6)

# MATLAB front-end
add_custom_target(matlab DEPENDS casadiMEX)

# Find MATLAB
find_package(MATLAB REQUIRED)

include_directories(${MATLAB_INCLUDE_DIR})
add_definitions(${MATLAB_EXTRA_CXXFLAGS})

# Generate SWIG wrapper
set_source_files_properties(../casadi.i  PROPERTIES  CPLUSPLUS ON)
swig_add_module(casadiMEX matlab ../casadi.i)
swig_link_libraries(casadiMEX ${MATLAB_LIBRARIES} casadi)
set_target_properties(casadiMEX PROPERTIES PREFIX "" SUFFIX .${MATLAB_MEX_EXT})

# Install the generated front-end to ${CMAKE_INSTALL_PREFIX}
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/+casadi DESTINATION ${CMAKE_INSTALL_PREFIX})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/SwigRef.m DESTINATION ${CMAKE_INSTALL_PREFIX})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/SwigMem.m DESTINATION ${CMAKE_INSTALL_PREFIX})
install(TARGETS casadiMEX DESTINATION ${CMAKE_INSTALL_PREFIX})

# Install matlab helper files
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/swig_typename_convertor_matlab2cpp.m
              ${CMAKE_CURRENT_SOURCE_DIR}/swig_typename_convertor_cpp2matlab.m
  DESTINATION ${CMAKE_INSTALL_PREFIX})
