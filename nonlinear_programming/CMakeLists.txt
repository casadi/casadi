cmake_minimum_required(VERSION 2.6)
SET(NONLINEAR_PROGRAMMING_SRCS
  symbolic_nlp.hpp         symbolic_nlp.cpp
  sqp_method.hpp          sqp_method.cpp          sqp_internal.hpp          sqp_internal.cpp
  stabilized_sqp_method.hpp          stabilized_sqp_method.cpp          stabilized_sqp_internal.hpp          stabilized_sqp_internal.cpp
  scpgen.hpp              scpgen.cpp              scpgen_internal.hpp       scpgen_internal.cpp
  nlp_qp_solver.hpp       nlp_qp_solver.cpp       nlp_qp_internal.hpp       nlp_qp_internal.cpp
  nlp_implicit_solver.hpp nlp_implicit_solver.cpp nlp_implicit_internal.hpp nlp_implicit_internal.cpp
  newton_implicit_solver.hpp newton_implicit_solver.cpp newton_implicit_internal.hpp newton_implicit_internal.cpp
  simple_homotopy_nlp_solver.hpp simple_homotopy_nlp_solver.cpp simple_homotopy_nlp_internal.hpp simple_homotopy_nlp_internal.cpp
)

if(ENABLE_STATIC)
add_library(casadi_nonlinear_programming STATIC ${NONLINEAR_PROGRAMMING_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(casadi_nonlinear_programming SHARED ${NONLINEAR_PROGRAMMING_SRCS})
target_link_libraries(casadi_nonlinear_programming casadi)
if(WITH_CSPARSE)
  target_link_libraries(casadi_nonlinear_programming casadi_csparse_interface)
endif(WITH_CSPARSE)

endif(ENABLE_SHARED)
install(TARGETS casadi_nonlinear_programming
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY ./ 
   DESTINATION include/casadi/nonlinear_programming
   FILES_MATCHING PATTERN "*.hpp"
   PATTERN ".svn" EXCLUDE
)

