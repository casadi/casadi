cmake_minimum_required(VERSION 2.0)

#find_package(SNOPT REQUIRED QUIET)

include_directories(${SNOPT_INCLUDE_DIRS})

SET(SNOPT_INTERFACE_SRCS
  SnoptInterface.cpp
  SnoptInterface.hpp
)

add_library(snopt_interface STATIC ${SNOPT_INTERFACE_SRCS})
if(ENABLE_SHARED)
add_library(snopt_interface_shared SHARED ${SNOPT_INTERFACE_SRCS})
endif(ENABLE_SHARED)
