cmake_minimum_required(VERSION 2.8.6)

# Generate SWIG wrapper (for both MATLAB and Octave)
set_source_files_properties(../casadi.i  PROPERTIES  CPLUSPLUS ON)
swig_module_initialize(casadiMEX matlab)
swig_add_source_to_module(casadiMEX FALSE MEX_FILE ../casadi.i)
get_directory_property(swig_extra_clean_files ADDITIONAL_MAKE_CLEAN_FILES)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${MEX_FILE}")

# Install the generated front-end to ${MATLAB_PREFIX}
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/+casadi DESTINATION ${MATLAB_PREFIX})

# Install MATLAB/Octave helper files
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/SwigRef.m DESTINATION ${MATLAB_PREFIX})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/SwigMem.m DESTINATION ${MATLAB_PREFIX})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/SwigGet.m DESTINATION ${MATLAB_PREFIX})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/SwigStorage.m DESTINATION ${MATLAB_PREFIX})
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/swig_typename_convertor_matlab2cpp.m
              ${CMAKE_CURRENT_SOURCE_DIR}/SwigType.m
              ${CMAKE_CURRENT_SOURCE_DIR}/returntypes_call.m
              ${CMAKE_CURRENT_SOURCE_DIR}/returntypes.m
  DESTINATION ${MATLAB_PREFIX})

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/OptiCallbackHelper.m
  DESTINATION ${MATLAB_PREFIX}/+casadi)
# MATLAB specific stuff
if (WITH_MATLAB)
  # Make target
  add_custom_target(matlab DEPENDS casadiMATLAB)

  # Find MATLAB
  find_package(MATLAB REQUIRED)

  # Compile MEX file
  if (EXISTS ${MATLAB_INCLUDE_DIR}/../version/cpp_mexapi_version.cpp)
    add_library(casadiMATLAB MODULE ${MEX_FILE} ${MATLAB_INCLUDE_DIR}/../version/cpp_mexapi_version.cpp)
  else()
    add_library(casadiMATLAB MODULE ${MEX_FILE})
  endif()
  target_include_directories(casadiMATLAB PUBLIC ${MATLAB_INCLUDE_DIR})
  target_link_libraries(casadiMATLAB ${MATLAB_LIBRARIES} casadi)
  set_target_properties(casadiMATLAB
    PROPERTIES OUTPUT_NAME "casadiMEX"
               PREFIX ""
               SUFFIX .${MATLAB_MEX_EXT})
  install(TARGETS casadiMATLAB DESTINATION ${MATLAB_PREFIX})
endif (WITH_MATLAB)

# Octave specific stuff
if (WITH_OCTAVE)
  # Make target
  add_custom_target(octave DEPENDS casadiOCTAVE)

  # Find Octave
  if (NOT OCTAVE_LIBRARIES)
  find_package(OCTAVE REQUIRED)
  endif()

  # Compile MEX file
  add_library(casadiOCTAVE MODULE ${MEX_FILE})
  target_include_directories(casadiOCTAVE PUBLIC ${OCTAVE_INCLUDE_DIRS})
  if(APPLE)
    target_link_libraries(casadiOCTAVE casadi)
    set_target_properties(casadiOCTAVE PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
  else()
    target_link_libraries(casadiOCTAVE ${OCTAVE_LIBRARIES} casadi)
  endif()
  set_target_properties(casadiOCTAVE
    PROPERTIES OUTPUT_NAME "casadiMEX"
               PREFIX ""
               SUFFIX .mex)
  install(TARGETS casadiOCTAVE DESTINATION ${MATLAB_PREFIX})
    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/../octave/import.m DESTINATION ${MATLAB_PREFIX}/octave-extra)
    if (WITH_OCTAVE_IMPORT)
      message(WARNING "CMake option WITH_OCTAVE_IMPORT has been deprecated. 'import.m' is now copied to the octave-extra subdirectory.")
      install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/../octave/import.m DESTINATION ${MATLAB_PREFIX})
  endif()
endif (WITH_OCTAVE)

# Example of how to extend CasADi with additional features
if (WITH_EXTENDING_CASADI AND WITH_MATLAB)
  set_source_files_properties(../extending_casadi/extending_casadi.i PROPERTIES CPLUSPLUS ON)
  swig_add_module(extending_casadiMEX matlab ../extending_casadi/extending_casadi.i)
  swig_link_libraries(extending_casadiMEX ${MATLAB_LIBRARIES} extending_casadi)
  set_target_properties(extending_casadiMEX PROPERTIES PREFIX "" SUFFIX .${MATLAB_MEX_EXT})

  # Install the generated front-end to ${MATLAB_PREFIX}
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/+extending_casadi DESTINATION ${MATLAB_PREFIX})
  install(TARGETS extending_casadiMEX DESTINATION ${MATLAB_PREFIX})

  add_custom_target(extending_casadi_matlab DEPENDS extending_casadiMEX extending_casadi)
endif()
