cmake_minimum_required(VERSION 2.6)
SET(CONVEX_PROGRAMMING_SRCS
 qp_lp_solver.cpp    qp_lp_solver.hpp    qp_lp_internal.cpp    qp_lp_internal.hpp
 sdp_socp_solver.cpp sdp_socp_solver.hpp sdp_socp_internal.cpp sdp_socp_internal.hpp
)

if(ENABLE_STATIC)
add_library(casadi_convex_programming STATIC ${CONVEX_PROGRAMMING_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(casadi_convex_programming SHARED ${CONVEX_PROGRAMMING_SRCS})
endif(ENABLE_SHARED)
install(TARGETS casadi_convex_programming
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY ./ 
   DESTINATION include/casadi/convex_programming
   FILES_MATCHING PATTERN "*.hpp"
   PATTERN ".svn" EXCLUDE
)

