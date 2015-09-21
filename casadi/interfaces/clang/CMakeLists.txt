cmake_minimum_required(VERSION 2.8.6)
include_directories(${CLANG_INCLUDE_DIR})
add_definitions(${CLANG_DEFINITIONS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CLANG_CXX_FLAGS}")

casadi_plugin(Compiler clang
  clang_compiler.hpp
  clang_compiler.cpp
  clang_compiler_meta.cpp)

casadi_plugin_link_libraries(Compiler clang ${CLANG_LIBRARIES})

if(MINGW)
  find_library(SHLWAPI libshlwapi /usr/${PREFIX}/lib)
  casadi_plugin_link_libraries(Compiler clang ${SHLWAPI})
endif()

# Get resource directory
file(COPY ${CLANG_LLVM_LIB_DIR}/clang DESTINATION ${LIBRARY_OUTPUT_PATH}/casadi/jit)

# Create text files containing system paths
file(WRITE system_includes.txt) # for both C and C++
file(WRITE csystem_includes.txt) # Only C
file(WRITE cxxsystem_includes.txt) # Only C++

# Create text files containing default options
file(WRITE system_flags.txt) # for both C and C++
file(WRITE csystem_flags.txt) # Only C
file(WRITE cxxsystem_flags.txt) # Only C++

# On Windows, use relative paths
if (WIN32)
  file(APPEND system_includes.txt "mingw\n")
  file(APPEND cxxsystem_includes.txt "mingw/c++\n")
else()
  # Get C and C++ system includes
  # Note on OSX, "gcc" is often an alias for clang
  set(GET_SYSTEM_INCLUDES_SH ${CMAKE_CURRENT_SOURCE_DIR}/../../../misc/get_system_includes.sh)
  execute_process(COMMAND ${GET_SYSTEM_INCLUDES_SH} gcc c
                  OUTPUT_VARIABLE CSYSTEM_INCLUDES)
  file(APPEND csystem_includes.txt "${CSYSTEM_INCLUDES}")
  execute_process(COMMAND ${GET_SYSTEM_INCLUDES_SH} gcc c++
                  OUTPUT_VARIABLE CXXSYSTEM_INCLUDES)
  file(APPEND cxxsystem_includes.txt "${CXXSYSTEM_INCLUDES}")

  # Printout
  message(STATUS "Clang JIT uses the following C includes:\n"
                 ${CSYSTEM_INCLUDES}
                 "  and the following C++ includes:\n"
                 ${CXXSYSTEM_INCLUDES}
                 "  To change this, manually edit csystem_paths.txt and cxxsystem_paths.txt\n"
                 "  located in the <installation dir>/lib/casadi/jit/ directory\n"
                 "  or use the accompanied get_system_includes.sh shell script.")

  # Copy the shell script to the same directory
  file(COPY ${GET_SYSTEM_INCLUDES_SH} DESTINATION ${LIBRARY_OUTPUT_PATH}/casadi/jit)
endif()

# Copy text files with paths to the designated location
file(COPY system_includes.txt csystem_includes.txt cxxsystem_includes.txt
          system_flags.txt csystem_flags.txt cxxsystem_flags.txt
     DESTINATION ${LIBRARY_OUTPUT_PATH}/casadi/jit)

if(WITH_SELFCONTAINED)
if (WITH_MATLAB)
install(DIRECTORY ${LIBRARY_OUTPUT_PATH}/casadi/jit
        DESTINATION ${CMAKE_INSTALL_PREFIX}/casadi)
else()
install(DIRECTORY ${LIBRARY_OUTPUT_PATH}/casadi/jit
        DESTINATION ${CMAKE_INSTALL_PREFIX}/casadi/casadi)
endif()
else()
install(DIRECTORY ${LIBRARY_OUTPUT_PATH}/casadi/jit
        DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/casadi)
endif()
