cmake_minimum_required(VERSION 2.6)

include_directories(include src/bounds src/lp src/sdp src/solver src/sys src/vecmat)

# cmake_policy use MACOSX_RPATH=ON by default
if(COMMAND cmake_policy)
  if(CMAKE_MAJOR_VERSION GREATER 2)
    cmake_policy(SET CMP0042 NEW)
  endif()
endif()

file(GLOB DSDP_SRC "src/bounds/*.c"  "src/bounds/*.h" "src/lp/*.c" "src/lp/*.h" "src/sdp/*.c" "src/sdp/*.h" "src/solver/*.c" "src/solver/*.h" "src/sys/*.c" "src/sys/*.h" "src/vecmat/*.c" "src/vecmat/*.h")

set(DSDP_LIBRARIES ${DSDP_LIBRARIES})

if(ENABLE_STATIC)
add_library(casadi_dsdp STATIC ${DSDP_SRC})
endif()
if(ENABLE_SHARED)
add_library(casadi_dsdp SHARED ${DSDP_SRC})
endif()
target_link_libraries(casadi_dsdp ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
install(TARGETS casadi_dsdp
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY include
  DESTINATION include/casadi/external_packages/DSDP
  FILES_MATCHING PATTERN "*.h"
  PATTERN ".svn" EXCLUDE)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
set_target_properties(
  casadi_dsdp
PROPERTIES
  SOVERSION ${PACKAGE_SO_VERSION}
)
endif()
