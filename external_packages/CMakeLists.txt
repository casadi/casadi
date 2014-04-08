cmake_minimum_required(VERSION 2.6)

macro(casadi_external_library name)
if(ENABLE_STATIC)
add_library(${name} STATIC ${ARGN})
endif()
if(ENABLE_SHARED)
add_compiler_export_flags()
add_library(${name} SHARED ${ARGN})
generate_export_header(${name})
endif()
set(CASADI_MODULES ${CASADI_MODULES} ${name} PARENT_SCOPE)
install(TARGETS ${name}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
endmacro()

add_subdirectory(tinyxml)

if(WITH_CSPARSE)
  add_subdirectory(CSparse)
endif()

if(WITH_SUNDIALS)
  add_subdirectory(sundials-2.5mod)
endif()

if(QPOASES_FOUND)
  add_subdirectory(qpOASES)
endif()

if(DSDP_FOUND)
  add_subdirectory(DSDP)
endif()
