cmake_minimum_required(VERSION 3.10.2)

if(SWIG_EXPORT)
  set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_SOURCE_DIR}/target/extra)
  set(CMAKE_SWIG_WRAP_OUTDIR ${CMAKE_CURRENT_SOURCE_DIR}/target/source/dummy)
endif()

# Main header files
if(NOT SKIP_CONFIG_H_GENERATION)
  configure_file(../swig_config.h.cmake "${CMAKE_CURRENT_BINARY_DIR}/swig_config.h" ESCAPE_QUOTES)
endif()

if(SWIG_IMPORT)
  set(MEX_FILE ${CMAKE_CURRENT_SOURCE_DIR}/target/source/casadiMATLAB_wrap.cxx)
  set_source_files_properties(${MEX_FILE}  PROPERTIES  CPLUSPLUS ON)
else()
  # Generate SWIG wrapper (for both MATLAB and Octave)
  set_source_files_properties(../casadi.i  PROPERTIES  CPLUSPLUS ON)
  swig_module_initialize(casadiMEX matlab)
  swig_add_source_to_module(casadiMEX FALSE MEX_FILE ../casadi.i)
  get_directory_property(swig_extra_clean_files ADDITIONAL_MAKE_CLEAN_FILES)
  set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${MEX_FILE}")
endif()

if(SWIG_EXPORT)
  add_custom_target(matlab_source DEPENDS ${MEX_FILE})
endif()

set(SWIG_EXTRA ${CMAKE_CURRENT_BINARY_DIR})
if(SWIG_IMPORT)
  set(SWIG_EXTRA ${CMAKE_CURRENT_SOURCE_DIR}/target/extra)
endif()


# Install the generated front-end to ${MATLAB_PREFIX}
install(DIRECTORY ${SWIG_EXTRA}/+casadi DESTINATION ${MATLAB_PREFIX})

# Install MATLAB/Octave helper files
install(FILES ${SWIG_EXTRA}/SwigRef.m DESTINATION ${MATLAB_PREFIX})
install(FILES ${SWIG_EXTRA}/SwigMem.m DESTINATION ${MATLAB_PREFIX})
install(FILES ${SWIG_EXTRA}/SwigGet.m DESTINATION ${MATLAB_PREFIX})
install(FILES ${SWIG_EXTRA}/SwigStorage.m DESTINATION ${MATLAB_PREFIX})

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/swig_typename_convertor_matlab2cpp.m
              ${CMAKE_CURRENT_SOURCE_DIR}/SwigType.m
              ${CMAKE_CURRENT_SOURCE_DIR}/returntypes_call.m
              ${CMAKE_CURRENT_SOURCE_DIR}/returntypes.m
  DESTINATION ${MATLAB_PREFIX})

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/OptiCallbackHelper.m
  DESTINATION ${MATLAB_PREFIX}/+casadi)

# MATLAB specific stuff
if (WITH_MATLAB AND NOT SWIG_EXPORT)
  # Make target
  add_custom_target(matlab DEPENDS casadiMATLAB)

  # Find MATLAB
  find_package(MATLAB REQUIRED)

  # Compile MEX file
  add_library(casadiMATLAB MODULE ${MEX_FILE} ${MEX_VERSION_FILE})
  target_link_libraries(casadiMATLAB matlab::matlab casadi)
  if(MSVC)
    target_compile_options(casadiMATLAB PRIVATE "/bigobj")
  endif()
  set_target_properties(casadiMATLAB
    PROPERTIES OUTPUT_NAME "casadiMEX"
               PREFIX ""
               SUFFIX .${MATLAB_MEX_EXT})
  install(TARGETS casadiMATLAB DESTINATION ${MATLAB_PREFIX})
endif (WITH_MATLAB AND NOT SWIG_EXPORT)

# Octave specific stuff
if (WITH_OCTAVE AND NOT SWIG_EXPORT)
  # Make target
  add_custom_target(octave DEPENDS casadiOCTAVE)

  # Find Octave
  if (NOT OCTAVE_LIBRARIES)
  find_package(OCTAVE REQUIRED)
  endif()

  # Compile MEX file
  add_library(casadiOCTAVE MODULE ${MEX_FILE})
  target_compile_definitions(casadiOCTAVE PRIVATE HAVE_OCTAVE)
  target_link_libraries(casadiOCTAVE octave::octave casadi)
  if(APPLE)
    set_target_properties(casadiOCTAVE PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
  endif()
  if(MSVC)
    target_compile_options(casadiOCTAVE PRIVATE "/bigobj")
  endif()
  set_target_properties(casadiOCTAVE
    PROPERTIES OUTPUT_NAME "casadiMEX"
               PREFIX ""
               SUFFIX .${OCTAVE_MEX_EXT})

  install(TARGETS casadiOCTAVE DESTINATION ${MATLAB_PREFIX})
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/../octave/import.m DESTINATION ${MATLAB_PREFIX}/octave-extra)
  if (WITH_OCTAVE_IMPORT)
    message(WARNING "CMake option WITH_OCTAVE_IMPORT has been deprecated. 'import.m' is now copied to the octave-extra subdirectory.")
    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/../octave/import.m DESTINATION ${MATLAB_PREFIX})
  endif()
  if (TARGET octave::octinterp_adaptor)
    get_target_property(OCTINTERP_LOCATION octave::octinterp_adaptor IMPORTED_LOCATION_RELEASE)
    install(FILES ${OCTINTERP_LOCATION} DESTINATION ${LIB_PREFIX})
  endif()
endif (WITH_OCTAVE AND NOT SWIG_EXPORT)

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
