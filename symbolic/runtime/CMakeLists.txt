cmake_minimum_required(VERSION 2.6)

set(EMBEDDED_SRC "")
file(GLOB SRC *.hpp)

foreach(file ${SRC})
  get_filename_component(FILE_BASENAME ${file} NAME_WE)
  string(REGEX REPLACE ".hpp$" "_embedded.hpp" FILE_NEWNAME "${FILENAME}")
  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${FILE_BASENAME}_embedded.hpp"
    COMMAND ${CMAKE_COMMAND} -D "FILE_BASENAME=${FILE_BASENAME}" -D "FILE=${file}"  -D "BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}" -P ${CMAKE_CURRENT_SOURCE_DIR}/embedded.cmake 
    DEPENDS
      ${CMAKE_CURRENT_SOURCE_DIR}/embedded.cmake 
      ${CMAKE_CURRENT_SOURCE_DIR}/${FILE_BASENAME}.hpp
    VERBATIM
  )
  set(EMBEDDED_SRC ${EMBEDDED_SRC} "${CMAKE_CURRENT_BINARY_DIR}/${FILE_BASENAME}_embedded.hpp")
endforeach()

add_custom_target(casadi_runtime_embedded
              DEPENDS ${EMBEDDED_SRC})
              
set(RUNTIME_EMBEDDED_SRC "${EMBEDDED_SRC}" PARENT_SCOPE)



