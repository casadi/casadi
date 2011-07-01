cmake_minimum_required(VERSION 2.0)
project(casadi)

SET(CASADI_SRCS
  # MISC
  casadi_limits.hpp
  casadi_operators.hpp
  casadi_exception.hpp
  casadi_math.hpp
  printable_object.hpp        printable_object.cpp      # Base class enabling printing a Python-style "description" as well as a shorter "representation" of a class
  shared_object.hpp           shared_object.cpp         # This base class implements the reference counting (garbage collection) framework used in CasADi
  generic_type.hpp            generic_type.cpp          # Generic type used for options and for compatibility with dynamically typed languages like Python
  generic_type_internal.hpp   generic_type_internal.cpp # Internal class for the same
  options_functionality.cpp   options_functionality.hpp # Functionality for getting and setting options of a derived class
  pre_c99_support.hpp         pre_c99_support.cpp       # Some functions that ensures compatibility with C++ compilers not supporting the C99 standard
  stl_vector_tools.hpp                                  # Set of useful functions for the vector template class in STL

  # Template class Matrix<>, implements a sparse Matrix with row compressed storage, designed to work well with symbolic data types (SX)
  matrix/matrix.hpp                                     # The Matrix template class
  matrix/crs_sparsity.hpp     matrix/crs_sparsity.cpp   # reference counted sparsity pattern (several matrices may share the same sparsity pattern)
  matrix/slice.hpp            matrix/slice.cpp          # A "slice" in a language such as Python, A[start:stop:step], or Matlab, A(start:step:stop)
  matrix/submatrix.hpp                                  # A reference to a block of the matrix to allow operations such as A(:,3) = ...
  matrix/nonzeros.hpp                                   # A reference to a set of nonzeros of the matrix to allow operations such as A[3] = ...
  matrix/matrix_tools.hpp     matrix/matrix_tools.cpp   # Set of functions
  matrix/sparsity_tools.hpp   matrix/sparsity_tools.cpp # Set of functions for sparsity

  # Directed, acyclic graph representation with scalar expressions
  sx/sx.hpp                  sx/sx.cpp                  # Public, smart pointer class, 
  sx/sx_node.hpp             sx/sx_node.cpp             # Base class for all the nodes
  sx/symbolic_sx_node.hpp                               # A symbolic SX variable 
  sx/constant_sx_node.hpp                               # A constant SX node
  sx/binary_sx_node.cpp      sx/binary_sx_node.hpp      # A binary _or_ unary operation
  sx/sx_tools.cpp            sx/sx_tools.hpp            # Set of functions

  # More general graph representation with sparse matrix expressions and function evaluations
  mx/mx.cpp                  mx/mx.hpp                  # Public, smart pointer class, 
  mx/mx_node.cpp             mx/mx_node.hpp             # Base class for all the MX nodes
  mx/mx_constant.cpp         mx/mx_constant.hpp         # A constant MX node
  mx/symbolic_mx_node.cpp    mx/symbolic_mx_node.hpp    # A symbolic MX variable
  mx/unary_op.cpp            mx/unary_op.hpp            # Unary operation
  mx/scalar_matrix_op.cpp    mx/scalar_matrix_op.hpp    # Scalar-matrix elementwise binary operation
  mx/matrix_scalar_op.cpp    mx/matrix_scalar_op.hpp    # Matrix-scalar elementwise binary operation
  mx/matrix_matrix_op.cpp    mx/matrix_matrix_op.hpp    # Matrix-matrix elementwise binary operation
  mx/multiplication.cpp      mx/multiplication.hpp      # Matrix multiplication
  mx/evaluation.cpp          mx/evaluation.hpp          # FX function evaluation
  mx/norm.cpp                mx/norm.hpp                # 1-norm, 2-norm and infinity-norm
  mx/mapping.cpp             mx/mapping.hpp             # Maps non-zeros
  mx/densification.cpp       mx/densification.hpp       # Makes a matrix dense
  mx/if_else_node.cpp        mx/if_else_node.hpp        # An if-else branch
  mx/multiple_output.cpp     mx/multiple_output.hpp     # Base class for nodes with multiple outputs
  mx/jacobian_reference.cpp  mx/jacobian_reference.hpp  # A reference to a Jacobian block
  mx/mx_tools.cpp            mx/mx_tools.hpp
  # A dynamically created function with AD capabilities
  fx/fx.hpp                  fx/fx.cpp                  fx/fx_internal.hpp                  fx/fx_internal.cpp
  fx/x_function.hpp          fx/x_function.cpp          fx/x_function_internal.hpp          fx/x_function_internal.cpp
  fx/sx_function.hpp         fx/sx_function.cpp         fx/sx_function_internal.hpp         fx/sx_function_internal.cpp
  fx/mx_function.hpp         fx/mx_function.cpp         fx/mx_function_internal.hpp         fx/mx_function_internal.cpp
  fx/c_function.hpp          fx/c_function.cpp          fx/c_function_internal.hpp          fx/c_function_internal.cpp
  fx/external_function.hpp   fx/external_function.cpp   fx/external_function_internal.hpp   fx/external_function_internal.cpp
  fx/jacobian.hpp            fx/jacobian.cpp            fx/jacobian_internal.hpp            fx/jacobian_internal.cpp
  fx/linear_solver.hpp       fx/linear_solver.cpp       fx/linear_solver_internal.hpp       fx/linear_solver_internal.cpp
  fx/implicit_function.hpp   fx/implicit_function.cpp   fx/implicit_function_internal.hpp   fx/implicit_function_internal.cpp
  fx/integrator.hpp          fx/integrator.cpp          fx/integrator_internal.hpp          fx/integrator_internal.cpp
  fx/nlp_solver.hpp          fx/nlp_solver.cpp          fx/nlp_solver_internal.hpp          fx/nlp_solver_internal.cpp
  fx/simulator.hpp           fx/simulator.cpp           fx/simulator_internal.hpp           fx/simulator_internal.cpp
  fx/parallelizer.hpp        fx/parallelizer.cpp        fx/parallelizer_internal.hpp        fx/parallelizer_internal.cpp
  fx/ocp_solver.hpp          fx/ocp_solver.cpp          fx/ocp_solver_internal.hpp          fx/ocp_solver_internal.cpp
  fx/qp_solver.hpp           fx/qp_solver.cpp           fx/qp_solver_internal.hpp           fx/qp_solver_internal.cpp
  fx/function_io.hpp         fx/function_io.cpp
  fx/fx_tools.hpp            fx/fx_tools.cpp
)
if(ENABLE_STATIC)
add_library(casadi STATIC ${CASADI_SRCS})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(casadi SHARED ${CASADI_SRCS})
endif(ENABLE_SHARED)
install(TARGETS casadi
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
