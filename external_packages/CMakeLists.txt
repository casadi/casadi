cmake_minimum_required(VERSION 2.6)

macro(casadi_external_library name)
  # cmake_policy use MACOSX_RPATH=ON by default
  if(COMMAND cmake_policy)
    if(CMAKE_MAJOR_VERSION GREATER 2)
      cmake_policy(SET CMP0042 NEW)
    endif()
  endif()

  if(ENABLE_STATIC)
    add_library(${name} STATIC ${ARGN})
  endif()
  if(ENABLE_SHARED)
    add_library(${name} SHARED ${ARGN})
  endif()
  generate_export_header(${name})
  set(CASADI_MODULES ${CASADI_MODULES} ${name} PARENT_SCOPE)
  install(TARGETS ${name}
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
  )
  set_target_properties(
    ${name}
  PROPERTIES
    SOVERSION ${PACKAGE_SO_VERSION}
  )
endmacro()

if(WITH_BUILD_TINYXML)
  add_subdirectory(tinyxml-2.6.2)
endif()

if(WITH_BUILD_CSPARSE)
  add_subdirectory(CSparse)
endif()

if(WITH_BUILD_SUNDIALS)
  add_subdirectory(sundials-2.5mod)
endif()

if(QPOASES_FOUND)
  add_subdirectory(qpOASES)
endif()

if(WITH_BUILD_DSDP AND DSDP_FOUND)
  add_subdirectory(DSDP)
endif()
