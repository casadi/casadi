cmake_minimum_required(VERSION 2.8.6)

# Make the shared library prefix and suffix available in C++

add_subdirectory(runtime)

set(CASADI_PUBLIC
  # MISC
  casadi_limits.hpp
  casadi_types.hpp
  casadi_logger.hpp
  casadi_file.hpp
  casadi_interrupt.hpp
  exception.hpp
  calculus.hpp
  global_options.hpp
  casadi_meta.hpp
  printable_object.hpp        # Interface class enabling printing a Python-style "description" as well as a shorter "representation" of a class
  shared_object.hpp           # This base class implements the reference counting (garbage collection) framework used in CasADi
  weak_ref.hpp                # Provides weak reference functionality (non-owning smart pointers)
  generic_type.hpp            # Generic type used for options and for compatibility with dynamically typed languages like Python
  options.hpp                 # Functionality for passing options to a class
  std_vector_tools.hpp        # Set of useful functions for the vector template class in STL
  timing.hpp
  polynomial.hpp              # Helper class for differentiating and integrating simple polynomials

  # Template class Matrix<>, implements a sparse Matrix with col compressed storage, designed to work well with symbolic data types (SX)
  generic_expression.hpp      # Base class for SXElem MX and Matrix<>
  generic_matrix.hpp          # Base class for Matrix<> and MX
  matrix.hpp                  # The Matrix template class
  sparsity_interface.hpp
  sparsity.hpp                # reference counted sparsity pattern
  slice.hpp                   # A "slice" in a language such as Python, A[start:stop:step], or Matlab, A(start:step:stop)
  submatrix.hpp               # A reference to a block of the matrix to allow operations such as A(:,3) = ...
  nonzeros.hpp                # A reference to a set of nonzeros of the matrix to allow operations such as A[3] = ...

  # Directed, acyclic graph representation with scalar expressions
  sx/sx_elem.hpp

  # More general graph representation with sparse matrix expressions and function evaluations
  mx.hpp

  # A dynamically created function with AD capabilities
  function/function.hpp
  function/callback.hpp
  function/external.hpp
  function/jit.hpp
  function/linsol.hpp
  function/rootfinder.hpp
  function/integrator.hpp
  function/nlpsol.hpp
  function/conic.hpp
  function/dple.hpp
  function/interpolant.hpp
  function/expm.hpp
  function/code_generator.hpp
  function/switch.hpp
  function/bspline.hpp
  function/map.hpp
  function/importer.hpp

  # MISC useful stuff
  misc/integration_tools.hpp
  misc/nlp_builder.hpp
  misc/xml_node.hpp
  misc/xml_file.hpp
  misc/variable.hpp
  misc/dae_builder.hpp

  # User include class with the most essential includes
  core.hpp
)

set(CASADI_INTERNAL
# For a so/dll with default hidden visibility, we need to instantiate some templates
# Since this compilation takes some time, we put it up-front,
# So that a make -jn will not be blocked
casadi_instantiations.cpp

# MISC
casadi_logger.cpp
casadi_file.cpp
casadi_interrupt.cpp
global_options.cpp
${PROJECT_BINARY_DIR}/casadi_meta.cpp
shared_object.cpp
weak_ref.cpp
generic_type.cpp
generic_type_internal.hpp
options.cpp
std_vector_tools.cpp
timing.cpp
polynomial.cpp

# Template class Matrix<>, implements a sparse Matrix with col compressed storage, designed to work well with symbolic data types (SX)
matrix.cpp
sparse_storage.hpp   sparse_storage_impl.hpp
sparsity.cpp  sparsity_internal.hpp   sparsity_internal.cpp
slice.cpp

# Directed, acyclic graph representation with scalar expressions
sx/sx_elem.cpp             # Symbolic expression class (scalar-valued atomics)
sx/sx_node.hpp             sx/sx_node.cpp             # Base class for all the nodes
sx/symbolic_sx.hpp                                    # A symbolic SXElem variable
sx/constant_sx.hpp                                    # A constant SXElem node
sx/unary_sx.hpp                                       # A unary operation
sx/binary_sx.hpp                                      # A binary operation

# More general graph representation with sparse matrix expressions and function evaluations
mx.cpp                  # Symbolic expression class (matrix-valued atomics)
mx_node.hpp             mx_node.cpp             # Base class for all the MX nodes
constant_mx.hpp         constant_mx.cpp         # A constant MX node
symbolic_mx.hpp         symbolic_mx.cpp         # A symbolic MX variable
unary_mx.hpp            unary_mx.cpp            # Unary operation
binary_mx.hpp           binary_mx_impl.hpp      # Binary operation
multiplication.hpp      multiplication.cpp      # Matrix multiplication
einstein.hpp            einstein.cpp            # Einstein product
solve.hpp               solve_impl.hpp          # Solve linear system of equations
casadi_call.hpp         casadi_call.cpp         # Function call
casadi_find.hpp         casadi_find.cpp         # Find first nonzero
norm.hpp                norm.cpp                # 1-norm, 2-norm and infinity-norm
transpose.hpp           transpose.cpp           # Transpose
concat.hpp              concat.cpp              # Concatenation
split.hpp               split.cpp               # Split (inverse of concatenation)
reshape.hpp             reshape.cpp             # Reshape
subref.hpp              subref.cpp              # Submatrix reference
subassign.hpp           subassign.cpp           # Submatrix assignment
getnonzeros.hpp         getnonzeros.cpp         # Get the nonzeros of matrix
setnonzeros.hpp         setnonzeros_impl.hpp    # Assign the nonzeros of a matrix to another matrix
project.hpp             project.cpp             # Sparse copy to another sparsity pattern
determinant.hpp         determinant.cpp         # Determinant
inverse.hpp             inverse.cpp             # Matrix inverse
dot.hpp                 dot.cpp                 # Inner product
bilin.hpp               bilin.cpp               # Bilinear form
rank1.hpp               rank1.cpp               # Rank-1 update
multiple_output.hpp     multiple_output.cpp     # Base class for nodes with multiple outputs
assertion.hpp           assertion.cpp           # Assertion
monitor.hpp             monitor.cpp             # Monitor
repmat.hpp              repmat.cpp              # RepMat

# A dynamically created function with AD capabilities
function/function.cpp
function/function_internal.hpp   function/function_internal.cpp   # Function object class (internal API)
function/oracle_function.hpp     function/oracle_function.cpp     # Specialization of FunctionInternal to hold an oracle
function/callback.cpp            # Interface for user-defined function classes (public API)
function/callback_internal.cpp   function/callback_internal.hpp   # Interface for user-defined function classes (internal API)
function/plugin_interface.hpp                                     # Plugin interface for Function
function/factory.hpp                                              # Helper class for derivative function generation
function/x_function.hpp                                           # Base class for SXFunction and MXFunction
function/sx_function.hpp         function/sx_function.cpp
function/mx_function.hpp         function/mx_function.cpp
function/external.cpp
function/jit.cpp
function/linsol.cpp              function/linsol_internal.hpp  function/linsol_internal.cpp
function/rootfinder_impl.hpp     function/rootfinder.cpp
function/integrator_impl.hpp     function/integrator.cpp
function/nlpsol.hpp              function/nlpsol_impl.hpp        function/nlpsol.cpp
function/conic_impl.hpp          function/conic.cpp
function/dple_impl.hpp           function/dple.cpp
function/interpolant_impl.hpp    function/interpolant.cpp
function/expm_impl.hpp           function/expm.cpp
function/code_generator.cpp
function/switch.hpp              function/switch.cpp
function/bspline.hpp             function/bspline.cpp
function/map.hpp                 function/map.cpp
function/importer.cpp            function/importer_internal.hpp function/importer_internal.cpp

# MISC useful stuff
misc/integration_tools.cpp
misc/nlp_builder.cpp
misc/xml_node.cpp
misc/xml_file.cpp                misc/xml_file_internal.hpp                misc/xml_file_internal.cpp
misc/variable.cpp
misc/dae_builder.cpp

# runtime
runtime/shared.cpp
)

set(CASADI_SRCS ${CASADI_PUBLIC} ${CASADI_INTERNAL})


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
