cmake_minimum_required(VERSION 2.8.6)
include_directories(${CLANG_INCLUDE_DIR})

casadi_plugin(JitFunction clang
  clang_interface.hpp
  clang_interface.cpp
  clang_interface_meta.cpp)

casadi_plugin_link_libraries(JitFunction clang ${CLANG_LIBRARIES})

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/helper_finddep.c" "#include <math.h>")
# Find out the dependencies of math.h
set(CMAKE_C_FLAGS_BACKUP "${CMAKE_C_FLAGS}")
set(CMAKE_C_FLAGS "-O3 -MD -MF ${CMAKE_CURRENT_BINARY_DIR}/helper_finddep.d")
TRY_COMPILE(RESULT "${CMAKE_CURRENT_BINARY_DIR}/temp"
  "${CMAKE_CURRENT_BINARY_DIR}/helper_finddep.c"
)
message("${RESULT}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS_BACKUP}")

# Slurp the dependecy list from the exported file
file(READ "${CMAKE_CURRENT_BINARY_DIR}/helper_finddep.d" helper_finddep_deps)

string(REPLACE "\\" "" helper_finddep_deps "${helper_finddep_deps}")

string(REGEX MATCHALL "([^\ ]+\\.h)" FOO_DEPS "${helper_finddep_deps}")

foreach(D ${FOO_DEPS})
  string(REPLACE "include" ";" FOO_DEPS_DISSECT "${D}")
  list(APPEND FOO_DEPS_DISSECT_LIST ${FOO_DEPS_DISSECT})
  list(GET FOO_DEPS_DISSECT_LIST -1 TARGET)
  string(SUBSTRING ${TARGET} 1 -1 TARGET)
  string(REGEX REPLACE "(x86_64|i686)-linux-gnu/" "" TARGET "${TARGET}")
  string(REGEX REPLACE "c\\+\\+/[\\d\\.]+/[\\w\\d\\-_]+" "" TARGET "${TARGET}")
  GET_FILENAME_COMPONENT(TARGETDIR "${TARGET}" PATH)
  install(FILES "${D}"
    DESTINATION "${CMAKE_INSTALL_PREFIX}/include/casadi/jit/${TARGETDIR}"
  )
  message("${D} : ${TARGET} -> ${TARGETDIR}")
endforeach()

