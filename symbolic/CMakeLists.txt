cmake_minimum_required(VERSION 2.6)

include_directories(${CMAKE_BINARY_DIR})

add_subdirectory(runtime)

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
  generic_type_internal.hpp                             # Internal class for the same
  options_functionality.cpp   options_functionality.hpp # Functionality for getting and setting options of a derived class
  stl_vector_tools.hpp        stl_vector_tools.cpp      # Set of useful functions for the vector template class in STL
  profiling.cpp               profiling.hpp

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
  mx/mx.hpp                  mx/mx.cpp                  # Public, smart pointer class, 
  mx/mx_node.hpp             mx/mx_node.cpp             # Base class for all the MX nodes
  mx/constant_mx.hpp         mx/constant_mx.cpp         # A constant MX node
  mx/symbolic_mx.hpp         mx/symbolic_mx.cpp         # A symbolic MX variable
  mx/unary_mx.hpp            mx/unary_mx.cpp            # Unary operation
  mx/binary_mx.hpp           mx/binary_mx_impl.hpp      # Binary operation
  mx/multiplication.hpp      mx/multiplication_impl.hpp # Matrix multiplication
  mx/solve.hpp               mx/solve_impl.hpp          # Solve linear system of equations
  mx/call_fx.hpp             mx/call_fx.cpp             # FX function evaluation
  mx/norm.hpp                mx/norm.cpp                # 1-norm, 2-norm and infinity-norm
  mx/transpose.hpp           mx/transpose.cpp           # Transpose
  mx/vertcat.hpp             mx/vertcat.cpp             # Vertical concatenation
  mx/vertsplit.hpp           mx/vertsplit.cpp           # Vertical split
  mx/reshape.hpp             mx/reshape.cpp             # Reshape
  mx/subref.hpp              mx/subref.cpp              # Submatrix reference
  mx/subassign.hpp           mx/subassign.cpp           # Submatrix assignment
  mx/getnonzeros.hpp         mx/getnonzeros.cpp         # Get the nonzeros of matrix
  mx/setnonzeros.hpp         mx/setnonzeros_impl.hpp    # Assign the nonzeros of a matrix to another matrix
  mx/set_sparse.hpp          mx/set_sparse.cpp          # Sparse copy
  mx/determinant.hpp         mx/determinant.cpp         # Determinant
  mx/inverse.hpp             mx/inverse.cpp             # Matrix inverse
  mx/inner_prod.hpp          mx/inner_prod.cpp          # Inner product
  mx/multiple_output.hpp     mx/multiple_output.cpp     # Base class for nodes with multiple outputs
  mx/mx_tools.hpp            mx/mx_tools.cpp

  # A dynamically created function with AD capabilities
  fx/io_interface.hpp                                   # Input/output interface for FX and FXInternal
  fx/io_scheme.hpp           fx/io_scheme.cpp           fx/io_scheme_internal.hpp           fx/io_scheme_internal.cpp 
  fx/fx.hpp                  fx/fx.cpp                  fx/fx_internal.hpp                  fx/fx_internal.cpp
  fx/x_function_internal.hpp                            # Base class for SXFunction and MXFunction
  fx/sx_function.hpp         fx/sx_function.cpp         fx/sx_function_internal.hpp         fx/sx_function_internal.cpp
  fx/mx_function.hpp         fx/mx_function.cpp         fx/mx_function_internal.hpp         fx/mx_function_internal.cpp
  fx/c_function.hpp          fx/c_function.cpp          fx/c_function_internal.hpp          fx/c_function_internal.cpp
  fx/external_function.hpp   fx/external_function.cpp   fx/external_function_internal.hpp   fx/external_function_internal.cpp
  fx/derivative.hpp          fx/derivative.cpp          fx/derivative_internal.hpp          fx/derivative_internal.cpp
  fx/linear_solver.hpp       fx/linear_solver.cpp       fx/linear_solver_internal.hpp       fx/linear_solver_internal.cpp
  fx/symbolic_qr.hpp         fx/symbolic_qr.cpp         fx/symbolic_qr_internal.hpp         fx/symbolic_qr_internal.cpp
  fx/implicit_function.hpp   fx/implicit_function.cpp   fx/implicit_function_internal.hpp   fx/implicit_function_internal.cpp
  fx/integrator.hpp          fx/integrator.cpp          fx/integrator_internal.hpp          fx/integrator_internal.cpp
  fx/nlp_solver.hpp          fx/nlp_solver.cpp          fx/nlp_solver_internal.hpp          fx/nlp_solver_internal.cpp
  fx/simulator.hpp           fx/simulator.cpp           fx/simulator_internal.hpp           fx/simulator_internal.cpp
  fx/control_simulator.hpp   fx/control_simulator.cpp   fx/control_simulator_internal.hpp   fx/control_simulator_internal.cpp
  fx/parallelizer.hpp        fx/parallelizer.cpp        fx/parallelizer_internal.hpp        fx/parallelizer_internal.cpp
  fx/ocp_solver.hpp          fx/ocp_solver.cpp          fx/ocp_solver_internal.hpp          fx/ocp_solver_internal.cpp
  fx/qp_solver.hpp           fx/qp_solver.cpp           fx/qp_solver_internal.hpp           fx/qp_solver_internal.cpp
  fx/sdp_solver.hpp          fx/sdp_solver.cpp          fx/sdp_solver_internal.hpp          fx/sdp_solver_internal.cpp
  fx/sdqp_solver.hpp         fx/sdqp_solver.cpp         fx/sdqp_solver_internal.hpp         fx/sdqp_solver_internal.cpp
  fx/socp_solver.hpp         fx/socp_solver.cpp         fx/socp_solver_internal.hpp         fx/socp_solver_internal.cpp
  fx/qcqp_solver.hpp         fx/qcqp_solver.cpp         fx/qcqp_solver_internal.hpp         fx/qcqp_solver_internal.cpp
  fx/lp_solver.hpp           fx/lp_solver.cpp           fx/lp_internal.hpp                  fx/lp_internal.cpp
  fx/fx_tools.hpp            fx/fx_tools.cpp
  fx/xfunction_tools.hpp     fx/xfunction_tools.cpp
  fx/code_generator.hpp      fx/code_generator.cpp

  # User include class with the most essential includes
  casadi.hpp
  
  fx/schemes_metadata.hpp fx/schemes_metadata.cpp fx/schemes_helpers.hpp
  
)

set_source_files_properties( ${RUNTIME_EMBEDDED_SRC} PROPERTIES GENERATED TRUE )

if(ENABLE_STATIC)
add_library(casadi STATIC ${CASADI_SRCS} ${RUNTIME_EMBEDDED_SRC})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(casadi SHARED ${CASADI_SRCS} ${RUNTIME_EMBEDDED_SRC})
endif(ENABLE_SHARED)

add_dependencies(casadi casadi_runtime_embedded)

install(TARGETS casadi
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY ./                                                                                                                                                                                          
   DESTINATION include/casadi/symbolic                                                                                                                                                                          
   FILES_MATCHING PATTERN "*.hpp"                                                                                                                                                                             
   PATTERN ".svn" EXCLUDE                                                                                                                                                                                     
)                                                                                                                                                                                                             
                 
