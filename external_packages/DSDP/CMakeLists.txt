cmake_minimum_required(VERSION 2.6)

include_directories(include src/bounds src/lp src/sdp src/solver src/sys src/vecmat)

FILE(GLOB DSDP_SRC "src/bounds/*.c"  "src/bounds/*.h" "src/lp/*.c" "src/lp/*.h" "src/sdp/*.c" "src/sdp/*.h" "src/solver/*.c" "src/solver/*.h" "src/sys/*.c" "src/sys/*.h" "src/vecmat/*.c" "src/vecmat/*.h" )

SET(DSDP_LIBRARIES ${DSDP_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )

if(ENABLE_STATIC)
add_library(casadi_dsdp STATIC ${DSDP_SRC})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(casadi_dsdp SHARED ${DSDP_SRC})
endif(ENABLE_SHARED)
install(TARGETS casadi_dsdp
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY include 
   DESTINATION include/casadi/external_packages/DSDP
   FILES_MATCHING PATTERN "*.h"
   PATTERN ".svn" EXCLUDE
)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
set_target_properties(casadi_dsdp PROPERTIES COMPILE_FLAGS "-Wno-unused-value")
endif()
