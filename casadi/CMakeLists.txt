cmake_minimum_required(VERSION 2.0)
project(casadi)

#SET(CMAKE_CXX_FLAGS "-O3 -std=c++0x -DWITH_CPP0X")
# SET(CMAKE_CXX_FLAGS "-O3")

SET(CASADI_SRCS
  # MISC
  matrix_size.hpp
  casadi_exception.hpp
  printable_object.hpp
  printable_object.cpp
  shared_object.hpp
  shared_object.cpp
  option.hpp
  option.cpp
  options_functionality.cpp
  options_functionality.hpp
  element.hpp
  matrix.hpp
  crs_sparsity.hpp
  crs_sparsity.cpp

  sx/sx.hpp
  sx/sx.cpp
  sx/sx_node.hpp
  sx/sx_node.cpp
  sx/symbolic_sx_node.hpp
  sx/constant_sx_node.hpp
  sx/binary_sx_node.cpp
  sx/binary_sx_node.hpp
  sx/sx_matrix.hpp
  sx/sx_matrix.cpp

  mx/mx.cpp
  mx/mx.hpp
  mx/mx_node.cpp
  mx/mx_node.hpp
  mx/norm.cpp
  mx/norm.hpp
  mx/transpose.cpp
  mx/transpose.hpp
  mx/mx_constant.cpp
  mx/mx_constant.hpp
  mx/evaluation.cpp
  mx/evaluation.hpp
  mx/if_else_node.cpp
  mx/if_else_node.hpp
  mx/unary_op.cpp
  mx/unary_op.hpp
  mx/scalar_matrix_op.cpp
  mx/scalar_matrix_op.hpp
  mx/matrix_scalar_op.cpp
  mx/matrix_scalar_op.hpp
  mx/matrix_matrix_op.cpp
  mx/matrix_matrix_op.hpp
  mx/multiplication.cpp
  mx/multiplication.hpp
  mx/matrix_element.cpp
  mx/matrix_element.hpp
  mx/symbolic_mx_node.cpp
  mx/symbolic_mx_node.hpp  
  mx/vertcat.cpp
  mx/vertcat.hpp
  mx/horzcat.cpp
  mx/horzcat.hpp

  fx/fx.hpp
  fx/fx.cpp
  fx/function_io.hpp
  fx/function_io.cpp
  fx/c_function.hpp
  fx/c_function.cpp
  fx/external_function.hpp
  fx/external_function.cpp
  fx/jacobian.hpp
  fx/jacobian.cpp
  fx/sx_function.hpp
  fx/sx_function.cpp
  fx/mx_function.hpp
  fx/mx_function.cpp
  fx/linear_solver.hpp 
  fx/linear_solver.cpp
  fx/integrator.hpp 
  fx/integrator.cpp
  fx/integrator_internal.hpp 
  fx/integrator_internal.cpp
  fx/nlp_solver.hpp 
  fx/nlp_solver.cpp
  fx/nlp_solver_internal.hpp 
  fx/nlp_solver_internal.cpp
  fx/simulator.hpp 
  fx/simulator.cpp
  fx/simulator_internal.hpp 
  fx/simulator_internal.cpp

  expression_tools.cpp
  expression_tools.hpp
  stl_vector_tools.hpp
)

add_library(casadi STATIC ${CASADI_SRCS})
if(ENABLE_SHARED)
add_library(casadi_shared SHARED ${CASADI_SRCS})
endif(ENABLE_SHARED)
