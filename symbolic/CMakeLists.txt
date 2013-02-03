cmake_minimum_required(VERSION 2.0)

SET(CASADI_SRCS
  # MISC
  casadi_limits.hpp
  casadi_types.hpp
  casadi_exception.hpp
  casadi_calculus.hpp
  casadi_math.hpp
  casadi_options.hpp          casadi_options.cpp
  casadi_meta.hpp             casadi_meta.cpp
  printable_object.hpp        printable_object.cpp      # Base class enabling printing a Python-style "description" as well as a shorter "representation" of a class
  shared_object.hpp           shared_object.cpp         # This base class implements the reference counting (garbage collection) framework used in CasADi
  weak_ref.hpp                weak_ref.cpp              # Provides weak reference functionality (non-owning smart pointers)
  generic_type.hpp            generic_type.cpp          # Generic type used for options and for compatibility with dynamically typed languages like Python
  generic_type_internal.hpp   generic_type_internal.cpp # Internal class for the same
  options_functionality.cpp   options_functionality.hpp # Functionality for getting and setting options of a derived class
  stl_vector_tools.hpp        stl_vector_tools.cpp      # Set of useful functions for the vector template class in STL

  # Template class Matrix<>, implements a sparse Matrix with row compressed storage, designed to work well with symbolic data types (SX)
  matrix/generic_expression.hpp                         # Base class for SX MX and Matrix<>
  matrix/generic_expression_tools.hpp                   # Tools for SX MX and Matrix<>
  matrix/generic_matrix.hpp                             # Base class for Matrix<> and MX
  matrix/generic_matrix_tools.hpp                       # Tools for both Matrix<> and MX
  matrix/matrix.hpp           matrix/matrix_impl.hpp    # The Matrix template class
  matrix/crs_sparsity.hpp     matrix/crs_sparsity.cpp   matrix/crs_sparsity_internal.hpp   matrix/crs_sparsity_internal.cpp    # reference counted sparsity pattern
  matrix/slice.hpp            matrix/slice.cpp          # A "slice" in a language such as Python, A[start:stop:step], or Matlab, A(start:step:stop)
  matrix/submatrix.hpp                                  # A reference to a block of the matrix to allow operations such as A(:,3) = ...
  matrix/nonzeros.hpp                                   # A reference to a set of nonzeros of the matrix to allow operations such as A[3] = ...
  matrix/matrix_tools.hpp     matrix/matrix_tools.cpp   # Set of functions
  matrix/sparsity_tools.hpp   matrix/sparsity_tools.cpp # Set of functions for sparsity

  # Directed, acyclic graph representation with scalar expressions
  sx/sx.hpp                  sx/sx.cpp                  # Public, smart pointer class, 
  sx/sx_node.hpp             sx/sx_node.cpp             # Base class for all the nodes
  sx/symbolic_sx.hpp                                    # A symbolic SX variable 
  sx/constant_sx.hpp                                    # A constant SX node
  sx/unary_sx.hpp                                       # A unary operation
  sx/binary_sx.hpp                                      # A binary operation
  sx/sx_tools.cpp            sx/sx_tools.hpp            # Set of functions

  # More general graph representation with sparse matrix expressions and function evaluations
  mx/mx.cpp                  mx/mx.hpp                  # Public, smart pointer class, 
  mx/mx_node.cpp             mx/mx_node.hpp             # Base class for all the MX nodes
  mx/constant_mx.cpp         mx/constant_mx.hpp         # A constant MX node
  mx/symbolic_mx.cpp         mx/symbolic_mx.hpp         # A symbolic MX variable
  mx/unary_mx.cpp            mx/unary_mx.hpp            # Unary operation
  mx/binary_mx.cpp           mx/binary_mx.hpp           # Binary operation
  mx/multiplication.cpp      mx/multiplication.hpp      # Matrix multiplication
  mx/solve.cpp               mx/solve.hpp               # Solve linear system of equations
  mx/evaluation_mx.cpp       mx/evaluation_mx.hpp       # FX function evaluation
  mx/norm.cpp                mx/norm.hpp                # 1-norm, 2-norm and infinity-norm
  mx/mapping.cpp             mx/mapping.hpp             # Maps non-zeros
  mx/densification.cpp       mx/densification.hpp       # Makes a matrix dense
  mx/multiple_output.cpp     mx/multiple_output.hpp     # Base class for nodes with multiple outputs
  mx/mx_tools.cpp            mx/mx_tools.hpp
  # A dynamically created function with AD capabilities
  fx/fx.hpp                  fx/fx.cpp                  fx/fx_internal.hpp                  fx/fx_internal.cpp
  fx/x_function_internal.hpp                            # Base class for SXFunction and MXFunction
  fx/sx_function.hpp         fx/sx_function.cpp         fx/sx_function_internal.hpp         fx/sx_function_internal.cpp
  fx/mx_function.hpp         fx/mx_function.cpp         fx/mx_function_internal.hpp         fx/mx_function_internal.cpp
  fx/c_function.hpp          fx/c_function.cpp          fx/c_function_internal.hpp          fx/c_function_internal.cpp
  fx/external_function.hpp   fx/external_function.cpp   fx/external_function_internal.hpp   fx/external_function_internal.cpp
  fx/derivative.hpp          fx/derivative.cpp          fx/derivative_internal.hpp          fx/derivative_internal.cpp
  fx/linear_solver.hpp       fx/linear_solver.cpp       fx/linear_solver_internal.hpp       fx/linear_solver_internal.cpp
  fx/implicit_function.hpp   fx/implicit_function.cpp   fx/implicit_function_internal.hpp   fx/implicit_function_internal.cpp
  fx/integrator.hpp          fx/integrator.cpp          fx/integrator_internal.hpp          fx/integrator_internal.cpp
  fx/nlp_solver.hpp          fx/nlp_solver.cpp          fx/nlp_solver_internal.hpp          fx/nlp_solver_internal.cpp
  fx/simulator.hpp           fx/simulator.cpp           fx/simulator_internal.hpp           fx/simulator_internal.cpp
  fx/control_simulator.hpp   fx/control_simulator.cpp   fx/control_simulator_internal.hpp   fx/control_simulator_internal.cpp
  fx/parallelizer.hpp        fx/parallelizer.cpp        fx/parallelizer_internal.hpp        fx/parallelizer_internal.cpp
  fx/ocp_solver.hpp          fx/ocp_solver.cpp          fx/ocp_solver_internal.hpp          fx/ocp_solver_internal.cpp
  fx/qp_solver.hpp           fx/qp_solver.cpp           fx/qp_solver_internal.hpp           fx/qp_solver_internal.cpp
  fx/sdp_solver.hpp          fx/sdp_solver.cpp          fx/sdp_solver_internal.hpp          fx/sdp_solver_internal.cpp
  fx/fx_tools.hpp            fx/fx_tools.cpp
  fx/xfunction_tools.hpp     fx/xfunction_tools.cpp

  # User include class with the most essential includes
  casadi.hpp
  
  autogenerated.hpp autogenerated.cpp

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
install(DIRECTORY ./                                                                                                                                                                                          
   DESTINATION include/casadi/symbolic                                                                                                                                                                          
   FILES_MATCHING PATTERN "*.hpp"                                                                                                                                                                             
   PATTERN ".svn" EXCLUDE                                                                                                                                                                                     
)                                                                                                                                                                                                             
                 
