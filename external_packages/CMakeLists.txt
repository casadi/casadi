cmake_minimum_required(VERSION 2.0)

add_subdirectory(tinyxml)
add_subdirectory(CSparse)
if(WITH_SUPERLU)
add_subdirectory(superlu_4_1)
endif(WITH_SUPERLU)
