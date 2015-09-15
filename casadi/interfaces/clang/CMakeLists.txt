cmake_minimum_required(VERSION 2.8.6)
include_directories(${CLANG_INCLUDE_DIR})
add_definitions(${LLVM_DEFINITIONS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CLANG_CXX_FLAGS}")

casadi_plugin(JitCompiler clang
  clang_interface.hpp
  clang_interface.cpp
  clang_interface_meta.cpp)

casadi_plugin_link_libraries(JitCompiler clang ${CLANG_LIBRARIES})

# Get C library
set(CLANG_CLIB)
list(APPEND CLANG_CLIB "#include <math.h>\n")

# Get C++ headers and the C headers they depend on
file(GLOB CLANG_CPPLIB ${LLVM_INSTALL_PREFIX}/include/c++/v1/*)
foreach(ITEM ${CLANG_CPPLIB})
   if( IS_DIRECTORY "${ITEM}" )
      INSTALL( DIRECTORY ${ITEM} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/casadi/jit/)
   else()
      INSTALL( FILES ${ITEM} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/casadi/jit/)
      FILE(READ "${ITEM}" LINES)
      foreach(L ${LINES})
        # Check if line is of form "#include <*.h>\n" and if so append to CLANG_CLIB
        string(REGEX MATCHALL "#include <[a-z_]+\\.h>\n" listinclude "${L}")
        foreach(LI ${listinclude})
          list(APPEND CLANG_CLIB "${LI}")
        endforeach()
      endforeach()
   endif()
endforeach()

# Remove duplicate includes and output to file
list (REMOVE_DUPLICATES CLANG_CLIB)

# Check if include is available and locate dependencies
set(CMAKE_C_FLAGS_BACKUP "${CMAKE_C_FLAGS}")
set(CMAKE_C_FLAGS "-O3 -MD -MF ${CMAKE_CURRENT_BINARY_DIR}/helper_finddep.d")
set(FOO_DEPS)
foreach(ITEM ${CLANG_CLIB})
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/helper_finddep.c" ${ITEM})
  try_compile(RESULT "${CMAKE_CURRENT_BINARY_DIR}/temp"
    "${CMAKE_CURRENT_BINARY_DIR}/helper_finddep.c")
  if (EXISTS ${CMAKE_CURRENT_BINARY_DIR}/helper_finddep.d)
    # Slurp the dependecy list from the exported file
    file(READ "${CMAKE_CURRENT_BINARY_DIR}/helper_finddep.d" helper_finddep_deps)
    string(REPLACE "\\" "" helper_finddep_deps "${helper_finddep_deps}")
    string(REGEX MATCHALL "([^\ ]+\\.h)" FOO_DEPS_ITEM "${helper_finddep_deps}")
    list(APPEND FOO_DEPS ${FOO_DEPS_ITEM})
  endif()
endforeach()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS_BACKUP}")

# Copy corresponding files
list (REMOVE_DUPLICATES FOO_DEPS)
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
