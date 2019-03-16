cmake_minimum_required(VERSION 2.8.6)

#set(EMBEDDED 1)
set(ENABLE_MKL_PARDISO OFF CACHE BOOL "Enable MKL Pardiso solver")
set(PROFILING OFF CACHE BOOL "Enable solver profiling (timing)")
set(CMAKE_INSTALL_INCLUDEDIR ${INCLUDE_PREFIX})
set(CMAKE_INSTALL_LIBDIR ${LIB_PREFIX})
set(CMAKE_INSTALL_BINDIR ${BIN_PREFIX})
add_subdirectory(osqp)

