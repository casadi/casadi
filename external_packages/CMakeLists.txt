cmake_minimum_required(VERSION 2.0)

add_subdirectory(tinyxml)
add_subdirectory(CSparse)

if(WITH_SUNDIALS)
  add_subdirectory(sundials-2.5mod)
endif()

if(WITH_QPOASES)
  add_subdirectory(qpOASES)
endif()
