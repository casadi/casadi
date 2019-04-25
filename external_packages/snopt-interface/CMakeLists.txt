cmake_minimum_required(VERSION 2.8.6)

project(snopt_interface C)

file(GLOB SRC snopt-interface/src/*.c)
file(GLOB SRC_FORTRAN snopt-interface/src/*.f90)
if(WIN32)
  add_library(casadi_snopt_interface STATIC ${SRC})
else()
  add_library(casadi_snopt_interface STATIC ${SRC} ${SRC_FORTRAN})
endif()
target_link_libraries(casadi_snopt_interface ${SNOPT_LIBRARIES})

target_include_directories(casadi_snopt_interface
    PUBLIC 
        $<INSTALL_INTERFACE:snopt-interface/include>    
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/snopt-interface/include>
)

