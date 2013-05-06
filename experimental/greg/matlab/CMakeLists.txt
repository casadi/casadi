cmake_minimum_required(VERSION 2.6)

find_package(MATLAB REQUIRED)
include_directories(${MATLAB_INCLUDE_DIRS})

SET(MATLAB_SRCS
#  MatlabInterface.cpp
#  MatlabInterface.hpp
	 engdemo.cpp
	 engdemo.hpp
)

add_library(matlab STATIC ${MATLAB_SRCS})

target_link_libraries(matlab ${MATLAB_LIBRARIES})

# add_library(snopt_interface STATIC ${SNOPT_INTERFACE_SRCS})
if(ENABLE_SHARED)
add_library(matlab_shared SHARED ${MATLAB_SRCS})
endif(ENABLE_SHARED)
