cmake_minimum_required(VERSION 3.10.2)

# Find MATLAB
find_package(MATLAB REQUIRED)
  
add_library(matlab_ipc SHARED
  matlab_ipc_external.c
)

target_compile_features(matlab_ipc PUBLIC c_std_99)

target_link_libraries(matlab_ipc matlab::matlab)

create_import_library(matlab_ipc matlab_ipc)

install(TARGETS matlab_ipc DESTINATION ${LIB_PREFIX})
