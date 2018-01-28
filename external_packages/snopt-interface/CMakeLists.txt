cmake_minimum_required(VERSION 2.8.6)

project(snopt_interface C)

include_directories(include)
include_directories(mod)
file(GLOB SRC src/*.c)
file(GLOB SRC_FORTRAN src/*.f90)
add_library(casadi_snopt_interface STATIC ${SRC} ${SRC_FORTRAN})
target_link_libraries(casadi_snopt_interface ${SNOPT_LIBRARIES})
