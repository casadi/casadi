cmake_minimum_required(VERSION 2.8.6)

# Make the shared library prefix and suffix available in C++
add_definitions(-DSHARED_LIBRARY_PREFIX="${CMAKE_SHARED_LIBRARY_PREFIX}")
add_definitions(-DSHARED_LIBRARY_SUFFIX="${CMAKE_SHARED_LIBRARY_SUFFIX}")

add_subdirectory(runtime)

set(CASADI_SRCS

  # For a so/dll with default hidden visibility, we need to instantiate some templates
  # Since this compilation takes some time, we put it up-front,
  # So that a make -jn will not be blocked
  casadi_instantiations.cpp
  
  # MISC
  casadi_limits.hpp
  casadi_types.hpp
  casadi_logger.hpp           casadi_logger.cpp
  casadi_interrupt.hpp        casadi_interrupt.cpp
  casadi_exception.hpp
  casadi_calculus.hpp
  casadi_math.hpp
  casadi_options.hpp          casadi_options.cpp
  casadi_meta.hpp             ${PROJECT_BINARY_DIR}/casadi_meta.cpp
  printable_object.hpp                                  # Interface class enabling printing a Python-style "description" as well as a shorter "representation" of a class
  shared_object.hpp           shared_object.cpp         # This base class implements the reference counting (garbage collection) framework used in CasADi
  weak_ref.hpp                weak_ref.cpp              # Provides weak reference functionality (non-owning smart pointers)
  generic_type.hpp            generic_type.cpp          # Generic type used for options and for compatibility with dynamically typed languages like Python
  generic_type_internal.hpp                             # Internal class for the same
  options_functionality.cpp   options_functionality.hpp # Functionality for getting and setting options of a derived class
  std_vector_tools.hpp        std_vector_tools.cpp      # Set of useful functions for the vector template class in STL
  profiling.hpp               profiling.cpp
  timing.hpp                  timing.cpp
  polynomial.hpp              polynomial.cpp            # Helper class for differentiating and integrating simple polynomials

  # Template class Matrix<>, implements a sparse Matrix with col compressed storage, designed to work well with symbolic data types (SX)
  matrix/generic_expression.hpp                         # Base class for SXElem MX and Matrix<>
  matrix/generic_matrix.hpp                             # Base class for Matrix<> and MX
  matrix/matrix.hpp           matrix/matrix_impl.hpp   matrix/matrix.cpp    # The Matrix template class
  matrix/sparse_storage.hpp   matrix/sparse_storage_impl.hpp
  matrix/sparsity_interface.hpp
  matrix/sparsity.hpp         matrix/sparsity.cpp  matrix/sparsity_internal.hpp   matrix/sparsity_internal.cpp    # reference counted sparsity pattern
  matrix/slice.hpp            matrix/slice.cpp          # A "slice" in a language such as Python, A[start:stop:step], or Matlab, A(start:step:stop)
  matrix/submatrix.hpp                                  # A reference to a block of the matrix to allow operations such as A(:,3) = ...
  matrix/nonzeros.hpp                                   # A reference to a set of nonzeros of the matrix to allow operations such as A[3] = ...

  # Directed, acyclic graph representation with scalar expressions
  sx/sx_elem.hpp             sx/sx_elem.cpp             # Symbolic expression class (scalar-valued atomics)
  sx/sx_node.hpp             sx/sx_node.cpp             # Base class for all the nodes
  sx/symbolic_sx.hpp                                    # A symbolic SXElem variable
  sx/constant_sx.hpp                                    # A constant SXElem node
  sx/unary_sx.hpp                                       # A unary operation
  sx/binary_sx.hpp                                      # A binary operation

  # More general graph representation with sparse matrix expressions and function evaluations
  mx/mx.hpp                  mx/mx.cpp                  # Symbolic expression class (matrix-valued atomics)
  mx/mx_node.hpp             mx/mx_node.cpp             # Base class for all the MX nodes
  mx/constant_mx.hpp         mx/constant_mx.cpp         # A constant MX node
  mx/symbolic_mx.hpp         mx/symbolic_mx.cpp         # A symbolic MX variable
  mx/unary_mx.hpp            mx/unary_mx.cpp            # Unary operation
  mx/binary_mx.hpp           mx/binary_mx_impl.hpp      # Binary operation
  mx/multiplication.hpp      mx/multiplication.cpp      # Matrix multiplication
  mx/solve.hpp               mx/solve_impl.hpp          # Solve linear system of equations
  mx/casadi_call.hpp         mx/casadi_call.cpp         # Function call
  mx/casadi_find.hpp         mx/casadi_find.cpp         # Find first nonzero
  mx/norm.hpp                mx/norm.cpp                # 1-norm, 2-norm and infinity-norm
  mx/transpose.hpp           mx/transpose.cpp           # Transpose
  mx/concat.hpp              mx/concat.cpp              # Concatenation
  mx/split.hpp               mx/split.cpp               # Split (inverse of concatenation)
  mx/reshape.hpp             mx/reshape.cpp             # Reshape
  mx/subref.hpp              mx/subref.cpp              # Submatrix reference
  mx/subassign.hpp           mx/subassign.cpp           # Submatrix assignment
  mx/getnonzeros.hpp         mx/getnonzeros.cpp         # Get the nonzeros of matrix
  mx/setnonzeros.hpp         mx/setnonzeros_impl.hpp    # Assign the nonzeros of a matrix to another matrix
  mx/project.hpp             mx/project.cpp             # Sparse copy to another sparsity pattern
  mx/determinant.hpp         mx/determinant.cpp         # Determinant
  mx/inverse.hpp             mx/inverse.cpp             # Matrix inverse
  mx/inner_prod.hpp          mx/inner_prod.cpp          # Inner product
  mx/multiple_output.hpp     mx/multiple_output.cpp     # Base class for nodes with multiple outputs
  mx/assertion.hpp           mx/assertion.cpp           # Assertion
  mx/monitor.hpp             mx/monitor.cpp             # Monitor
  mx/repmat.hpp              mx/repmat.cpp              # RepMat

  # A dynamically created function with AD capabilities
  function/io_interface.hpp                                         # Input/output interface for Function and FunctionInternal
  function/function.hpp            function/function.cpp            function/function_internal.hpp            function/function_internal.cpp
  function/plugin_interface.hpp                                     # Plugin interface for Function
  function/x_function.hpp                                           # Base class for SXFunction and MXFunction
  function/sx_function.hpp         function/sx_function.cpp
  function/mx_function.hpp         function/mx_function.cpp
  function/external.hpp            function/external.cpp
  function/linsol.hpp              function/linsol.cpp
  function/rootfinder.hpp          function/rootfinder.cpp
  function/integrator.hpp          function/integrator.cpp
  function/nlpsol.hpp          function/nlpsol.cpp
  function/simulator.hpp           function/simulator.cpp           function/simulator_internal.hpp           function/simulator_internal.cpp
  function/qp_solver.hpp           function/qp_solver.cpp
  function/code_generator.hpp      function/code_generator.cpp
  function/wrapper.hpp             function/adaptor.hpp
  function/switch.hpp              function/switch.cpp
  function/map.hpp                 function/map.cpp
  function/mapaccum.hpp            function/mapaccum.cpp
  function/compiler.hpp            function/compiler.cpp            function/compiler_internal.hpp function/compiler_internal.cpp
  function/kernel_sum.hpp          function/kernel_sum.cpp
  function/callback.cpp            function/callback.hpp            function/callback_internal.cpp function/callback_internal.hpp

  # MISC useful stuff
  misc/integration_tools.hpp       misc/integration_tools.cpp
  misc/nlp_builder.hpp             misc/nlp_builder.cpp
  misc/xml_node.hpp                misc/xml_node.cpp
  misc/xml_file.hpp                misc/xml_file.cpp                misc/xml_file_internal.hpp                misc/xml_file_internal.cpp
  misc/variable.hpp                misc/variable.cpp
  misc/dae_builder.hpp             misc/dae_builder.cpp

  # User include class with the most essential includes
  core.hpp
)

set_source_files_properties( ${RUNTIME_EMBEDDED_SRC} PROPERTIES GENERATED TRUE )

casadi_library(casadi ${CASADI_SRCS} ${RUNTIME_EMBEDDED_SRC})

add_dependencies(casadi casadi_runtime_embedded)

if(WITH_DL)
  # Core needs support for dynamic linking
  target_link_libraries(casadi ${CMAKE_DL_LIBS})
endif()

if(WITH_OPENCL)
  # Core depends on OpenCL for GPU calculations
  target_link_libraries(casadi ${OPENCL_LIBRARIES})
endif()

if(RT)
  # Realtime library
  target_link_libraries(casadi ${RT})
endif()

install(DIRECTORY ./
  DESTINATION include/casadi/core
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)
