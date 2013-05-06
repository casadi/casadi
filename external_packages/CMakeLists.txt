cmake_minimum_required(VERSION 2.6)

add_subdirectory(tinyxml)

if(WITH_CSPARSE)
  add_subdirectory(CSparse)
endif()

if(WITH_SUNDIALS)
  add_subdirectory(sundials-2.5mod)
endif()

if(WITH_QPOASES)
  add_subdirectory(qpOASES)
endif()

if(DSDP_FOUND)
  add_subdirectory(DSDP)
endif()
