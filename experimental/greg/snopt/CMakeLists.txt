cmake_minimum_required(VERSION 2.6)

#find_package(SNOPT REQUIRED QUIET)

include_directories(${SNOPT_INCLUDE_DIRS})

SET(SNOPT_INTERFACE_SRCS
  SnoptInterface.cpp
  SnoptInterface.hpp
)

if(ENABLE_STATIC)
add_library(snopt_interface STATIC ${SNOPT_INTERFACE_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(snopt_interface SHARED ${SNOPT_INTERFACE_SRCS})
endif(ENABLE_SHARED)
