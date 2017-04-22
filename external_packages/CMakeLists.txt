cmake_minimum_required(VERSION 2.8.6)

macro(casadi_external_library name)
  # cmake_policy use MACOSX_RPATH=ON by default
  if(COMMAND cmake_policy)
    if(CMAKE_MAJOR_VERSION GREATER 2)
      cmake_policy(SET CMP0042 NEW)
    endif()
  endif()
  # External packages are always build statically, not installed
  add_library(${name} STATIC ${ARGN})
endmacro()

if(WITH_TINYXML AND WITH_BUILD_TINYXML)
  add_subdirectory(tinyxml-2.6.2)
endif()

if(WITH_CSPARSE AND WITH_BUILD_CSPARSE)
  add_subdirectory(CSparse)
endif()

if(WITH_BUILD_SUNDIALS)
  add_subdirectory(casadi-sundials)
endif()

if(WITH_QPOASES)
  add_subdirectory(qpOASES)
endif()

if(WITH_BUILD_DSDP AND DSDP_FOUND)
  add_subdirectory(DSDP)
endif()

if(WITH_HPMPC AND WITH_BUILD_HPMPC)
  add_subdirectory(hpmpc)
endif()

if(WITH_BLASFEO AND WITH_BUILD_BLASFEO)
  add_subdirectory(blasfeo)
endif()
