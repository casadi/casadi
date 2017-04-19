cmake_minimum_required(VERSION 2.8.6)

# Turn a runtime file into a file with strings
macro(CASADI_STRINGIFY STRFILE)
  # Start with an empty file
  file(WRITE ${STRFILE} "")
  foreach(FILE ${ARGN})
    # Add declaration of string
    get_filename_component(FILENAME ${FILE} NAME_WE)
    file(APPEND ${STRFILE} "const char* ${FILENAME}_str =")
    # Append file as strings
    file(STRINGS ${FILE} FILE_CONTENTS)
    foreach(FILE_LINE ${FILE_CONTENTS})
      file(APPEND ${STRFILE} "\n  \"${FILE_LINE}\\n\"")
    endforeach()
    # End declaration
    file(APPEND ${STRFILE} ";\n\n")
  endforeach()
endmacro()

# Stringify C runtime
CASADI_STRINGIFY("${CMAKE_CURRENT_BINARY_DIR}/casadi_runtime_str.h"
  casadi_copy.hpp
  casadi_swap.hpp
  casadi_project.hpp
  casadi_densify.hpp
  casadi_sparsify.hpp
  casadi_scal.hpp
  casadi_axpy.hpp
  casadi_dot.hpp
  casadi_fill.hpp
  casadi_mtimes.hpp
  casadi_mv.hpp
  casadi_trans.hpp
  casadi_bilin.hpp)

set(EMBEDDED_SRC "")
file(GLOB SRC runtime.hpp)

foreach(file ${SRC})
  get_filename_component(FILE_BASENAME ${file} NAME_WE)
  string(REGEX REPLACE ".hpp$" "_embedded.hpp" FILE_NEWNAME "${FILENAME}")
  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${FILE_BASENAME}_embedded.hpp"
    COMMAND ${CMAKE_COMMAND} -D "FILE_BASENAME=${FILE_BASENAME}" -D "FILE=${file}"  -D "BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}" -P ${CMAKE_CURRENT_SOURCE_DIR}/embedded.cmake
    DEPENDS
      ${CMAKE_CURRENT_SOURCE_DIR}/embedded.cmake
      ${CMAKE_CURRENT_SOURCE_DIR}/${FILE_BASENAME}.hpp
    VERBATIM)
  set(EMBEDDED_SRC ${EMBEDDED_SRC} "${CMAKE_CURRENT_BINARY_DIR}/${FILE_BASENAME}_embedded.hpp")
endforeach()

add_custom_target(casadi_runtime_embedded DEPENDS ${EMBEDDED_SRC})

set(RUNTIME_EMBEDDED_SRC "${EMBEDDED_SRC}" PARENT_SCOPE)

install(FILES runtime.hpp shared.hpp
  DESTINATION include/casadi/core/runtime
)
