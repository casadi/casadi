cmake_minimum_required(VERSION 2.8.6)

#set(EMBEDDED 1)
set(ENABLE_MKL_PARDISO OFF CACHE BOOL "Enable MKL Pardiso solver")
set(PROFILING OFF CACHE BOOL "Enable solver profiling (timing)")
add_subdirectory(osqp)

