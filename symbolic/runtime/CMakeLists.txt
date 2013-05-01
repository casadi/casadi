cmake_minimum_required(VERSION 2.0)

set(EMBEDDED_SRC "")
file(GLOB SRC *.hpp)

foreach(file ${SRC})
  get_filename_component(FILE_BASENAME ${file} NAME_WE)
  string(REGEX REPLACE ".hpp$" "_embedded.hpp" FILE_NEWNAME "${FILENAME}")
  configure_file(embedded.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/${FILE_BASENAME}_embedded.cmake @ONLY)
  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${FILE_BASENAME}_embedded.hpp"
    COMMAND ${CMAKE_COMMAND} -D "FILE_BASENAME=${FILE_BASENAME}" -P ${CMAKE_CURRENT_BINARY_DIR}/${FILE_BASENAME}_embedded.cmake 
    DEPENDS
      ${CMAKE_CURRENT_BINARY_DIR}/${FILE_BASENAME}_embedded.cmake
      ${CMAKE_CURRENT_SOURCE_DIR}/${FILE_BASENAME}.hpp
    VERBATIM
  )
  set(EMBEDDED_SRC ${EMBEDDED_SRC} "${CMAKE_CURRENT_BINARY_DIR}/${FILE_BASENAME}_embedded.hpp")
endforeach()

if(ENABLE_STATIC)
add_library(casadi_cruntime_embedded STATIC ${EMBEDDED_SRC} dummy.cpp)
add_library(casadi_cruntime STATIC ${SRC} dummy.cpp)
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(casadi_cruntime_embedded SHARED ${EMBEDDED_SRC} dummy.cpp)
add_library(casadi_cruntime SHARED "runtime_h" ${SRC} dummy.cpp)
endif(ENABLE_SHARED)
install(TARGETS casadi_cruntime_embedded
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(TARGETS casadi_cruntime
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

install(DIRECTORY ./                                                                                                                                                                                          
   DESTINATION include/casadi/symbolic/cruntime                                                                                                                                                                          
   FILES_MATCHING PATTERN "*.h"                                                                                                                                                                             
   PATTERN ".svn" EXCLUDE                                                                                                                                                                                     
)
