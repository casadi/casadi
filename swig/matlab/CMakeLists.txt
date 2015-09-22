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
install(TARGETS casadiMEX DESTINATION ${CMAKE_INSTALL_PREFIX})

# Example of how to extend CasADi with additional features
if (WITH_EXTENDING_CASADI)
  set_source_files_properties(../extending_casadi/extending_casadi.i PROPERTIES CPLUSPLUS ON)
  swig_add_module(extending_casadiMEX matlab ../extending_casadi/extending_casadi.i)
  swig_link_libraries(extending_casadiMEX ${MATLAB_LIBRARIES} extending_casadi)
  set_target_properties(extending_casadiMEX PROPERTIES PREFIX "" SUFFIX .${MATLAB_MEX_EXT})

  # Install the generated front-end to ${CMAKE_INSTALL_PREFIX}
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/+extending_casadi DESTINATION ${CMAKE_INSTALL_PREFIX})
  install(TARGETS extending_casadiMEX DESTINATION ${CMAKE_INSTALL_PREFIX})

  add_custom_target(extending_casadi_matlab DEPENDS extending_casadiMEX extending_casadi)
endif()

# Install matlab helper files
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/SwigRef.m DESTINATION ${CMAKE_INSTALL_PREFIX})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/SwigMem.m DESTINATION ${CMAKE_INSTALL_PREFIX})
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/swig_typename_convertor_matlab2cpp.m
              ${CMAKE_CURRENT_SOURCE_DIR}/swig_typename_convertor_cpp2matlab.m
  DESTINATION ${CMAKE_INSTALL_PREFIX})
