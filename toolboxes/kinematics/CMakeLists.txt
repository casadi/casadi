cmake_minimum_required(VERSION 2.0)

include_directories(../..)

IF(ENABLE_STATIC)
   add_library(kinematics STATIC
     frame.hpp
     frame.cpp
     frame_node.hpp
     frame_node.cpp
     kinvec.hpp
     kinvec.cpp
     kinetics.hpp
     kinetics.cpp
   )
ENDIF(ENABLE_STATIC)

IF(ENABLE_SHARED)
   add_library(kinematics SHARED
     frame.hpp
     frame.cpp
     frame_node.hpp
     frame_node.cpp
     kinvec.hpp
     kinvec.cpp
     kinetics.hpp
     kinetics.cpp
   )
ENDIF(ENABLE_SHARED)

install(TARGETS kinematics
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
     

add_subdirectory(test)
add_subdirectory(examples)
