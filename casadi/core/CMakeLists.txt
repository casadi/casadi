cmake_minimum_required(VERSION 3.10.2)

# Create runtime
add_subdirectory(runtime)

# Binary interface to FMI, version 2.0
if(WITH_FMI2)
  # FMI C API, version 2
  add_definitions(-DWITH_FMI2)
  include_directories(${FMI2_INCLUDE_DIR})
  set(FMU2_SRC fmu2.hpp fmu2.cpp)
else()
  # Do not compile FMU 2 interface
  set(FMU2_SRC)
endif()

# Binary interface to FMI, version 3.0
if(WITH_FMI3)
  # FMI C API, version 2
  add_definitions(-DWITH_FMI3)
  include_directories(${FMI3_INCLUDE_DIR})
  set(FMU3_SRC fmu3.hpp fmu3.cpp)
else()
  # Do not compile FMU 2 interface
  set(FMU3_SRC)
endif()

set(CASADI_PUBLIC
  # MISC
  casadi_limits.hpp
  casadi_types.hpp
  casadi_common.hpp
  casadi_logger.hpp
  casadi_interrupt.hpp
  exception.hpp
  casadi_enum.hpp
  calculus.hpp
  global_options.hpp
  casadi_meta.hpp
  printable.hpp               # Interface class for printing to screen
  shared_object.hpp           # This base class implements the reference counting (garbage collection) framework used in CasADi
  generic_type.hpp            # Generic type used for options and for compatibility with dynamically typed languages like Python
  options.hpp                 # Functionality for passing options to a class
  casadi_misc.hpp             # Set of useful functions
  timing.hpp
  polynomial.hpp              # Helper class for differentiating and integrating simple polynomials

  # Template class Matrix<>, implements a sparse Matrix with col compressed storage, designed to work well with symbolic data types (SX)
  generic_expression.hpp      # Base class for SXElem MX and Matrix<>
  generic_matrix.hpp          # Base class for Matrix<> and MX
  generic_shared.hpp          # Base class for SharedObject
  generic_shared_internal.hpp # Base class for SharedObjectInternal
  generic_shared_impl.hpp
  matrix_fwd.hpp
  matrix_decl.hpp
  sx_fwd.hpp
  sx.hpp
  dm_fwd.hpp
  dm.hpp
  im_fwd.hpp
  im.hpp
  sparsity_interface.hpp
  sparsity.hpp                # reference counted sparsity pattern
  slice.hpp                   # A "slice" in a language such as Python, A[start:stop:step], or Matlab, A(start:step:stop)
  submatrix.hpp               # A reference to a block of the matrix to allow operations such as A(:,3) = ...
  nonzeros.hpp                # A reference to a set of nonzeros of the matrix to allow operations such as A[3] = ...

  # Directed, acyclic graph representation with scalar expressions
  sx_elem.hpp
  sx.hpp

  # More general graph representation with sparse matrix expressions and function evaluations
  mx.hpp

  # A dynamically created function with AD capabilities
  function.hpp
  callback.hpp
  external.hpp
  linsol.hpp
  rootfinder.hpp
  integrator.hpp
  nlpsol.hpp
  conic.hpp
  dple.hpp
  interpolant.hpp
  expm.hpp
  code_generator.hpp
  importer.hpp
  blazing_spline.hpp

  # MISC useful stuff
  integration_tools.hpp
  nlp_tools.hpp
  nlp_builder.hpp
  xml_node.hpp
  xml_file.hpp
  dae_builder.hpp
  optistack.hpp
  serializer.hpp
  serializing_stream.hpp
  fmu.hpp
  tools.hpp

  # User include class with the most essential includes
  core.hpp
)

set(CASADI_INTERNAL
  # MISC
  casadi_logger.cpp
  casadi_interrupt.cpp
  global_options.cpp
  ${CMAKE_CURRENT_BINARY_DIR}/../config.h
  casadi_meta.cpp
  generic_type.cpp
  generic_type_internal.hpp
  options.cpp
  casadi_misc.cpp
  casadi_common.cpp
  timing.cpp
  polynomial.cpp

  # Template class Matrix<>, implements a sparse Matrix with col compressed storage, designed to work well with symbolic data types (SX)
  matrix_impl.hpp
  sx_instantiator.cpp
  dm_instantiator.cpp
  im_instantiator.cpp

  sparsity.cpp  sparsity_internal.hpp   sparsity_internal.cpp
  slice.cpp generic_matrix.cpp
  shared_object.cpp

  # Directed, acyclic graph representation with scalar expressions
  sx_elem.cpp             # Symbolic expression class (scalar-valued atomics)
  sx_node.hpp             sx_node.cpp             # Base class for all the nodes
  symbolic_sx.hpp                                    # A symbolic SXElem variable
  constant_sx.hpp                                    # A constant SXElem node
  unary_sx.hpp                                       # A unary operation
  binary_sx.hpp                                      # A binary operation
  output_sx.hpp           output_sx.cpp
  call_sx.hpp

  # More general graph representation with sparse matrix expressions and function evaluations
  mx.cpp                  # Symbolic expression class (matrix-valued atomics)
  mx_node.hpp             mx_node.cpp             # Base class for all the MX nodes
  io_instruction.hpp      io_instruction.cpp      # Input and output MX instructions
  constant_mx.hpp         constant_mx.cpp         # A constant MX node
  symbolic_mx.hpp         symbolic_mx.cpp         # A symbolic MX variable
  unary_mx.hpp            unary_mx.cpp            # Unary operation
  binary_mx.hpp           binary_mx_impl.hpp      # Binary operation
  multiplication.hpp      multiplication.cpp      # Matrix multiplication
  einstein.hpp            einstein.cpp            # Einstein product
  solve.hpp               solve_impl.hpp          # Solve linear system of equations
  casadi_call.hpp         casadi_call.cpp         # Function call
  casadi_find.hpp         casadi_find.cpp         # Find first nonzero
  casadi_low.hpp          casadi_low.cpp          # Find first nonzero
  norm.hpp                norm.cpp                # 1-norm, 2-norm and infinity-norm
  mmin.hpp                mmin.cpp                # Matrix minimum
  transpose.hpp           transpose.cpp           # Transpose
  concat.hpp              concat.cpp              # Concatenation
  split.hpp               split.cpp               # Split (inverse of concatenation)
  reshape.hpp             reshape.cpp             # Reshape
  sparsity_cast.hpp       sparsity_cast.cpp       # Sparsity cast
  subref.hpp              subref.cpp              # Submatrix reference
  subassign.hpp           subassign.cpp           # Submatrix assignment
  getnonzeros.hpp         getnonzeros.cpp         # Get the nonzeros of matrix
  getnonzeros_param.hpp   getnonzeros_param.cpp   # Get the nonzeros of matrix, parametrically
  setnonzeros.hpp         setnonzeros_impl.hpp    # Assign the nonzeros of a matrix to another matrix
  setnonzeros_param.hpp   setnonzeros_param_impl.hpp  # Assign the nonzeros of a matrix to another matrix
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
  convexify.hpp           convexify.cpp           # Convexify
  logsumexp.hpp           logsumexp.cpp           # Logsumexp

  # A dynamically created function with AD capabilities
  function.cpp
  function_internal.hpp   function_internal.cpp   # Function object class (internal API)
  oracle_function.hpp     oracle_function.cpp     # Specialization of FunctionInternal to hold an oracle
  callback.cpp            # Interface for user-defined function classes (public API)
  callback_internal.cpp   callback_internal.hpp   # Interface for user-defined function classes (internal API)
  casadi_os.cpp           casadi_os.hpp           # Abstractions aroung operating system
  plugin_interface.hpp                                     # Plugin interface for Function
  factory.hpp                                              # Helper class for derivative function generation
  x_function.hpp                                           # Base class for SXFunction and MXFunction
  sx_function.hpp         sx_function.cpp
  mx_function.hpp         mx_function.cpp
  external_impl.hpp       external.cpp
  fmu_impl.hpp            fmu.cpp                  ${FMU2_SRC}  ${FMU3_SRC}
  fmu_function.hpp        fmu_function.cpp
  jit_function.hpp        jit_function.cpp
  linsol.cpp              linsol_internal.hpp  linsol_internal.cpp
  rootfinder_impl.hpp     rootfinder.cpp
  integrator_impl.hpp     integrator.cpp
  nlpsol.hpp              nlpsol_impl.hpp        nlpsol.cpp
  conic_impl.hpp          conic.cpp
  dple_impl.hpp           dple.cpp
  interpolant_impl.hpp    interpolant.cpp
  expm_impl.hpp           expm.cpp
  code_generator.cpp
  switch.hpp              switch.cpp
  bspline.hpp             bspline.cpp
  map.hpp                 map.cpp
  mapsum.hpp              mapsum.cpp
  finite_differences.hpp  finite_differences.cpp
  importer.cpp            importer_internal.hpp importer_internal.cpp
  blazing_spline.cpp blazing_spline_impl.hpp

  # MISC useful stuff
  integration_tools.cpp
  nlp_tools.cpp
  nlp_builder.cpp
  xml_node.cpp
  xml_file.cpp                xml_file_internal.hpp                xml_file_internal.cpp
  dae_builder.cpp             dae_builder_internal.hpp             dae_builder_internal.cpp
  optistack.cpp               optistack_internal.cpp               optistack_internal.hpp
  serializer.cpp              serializing_stream.cpp
  casadi_c.cpp
  tools.cpp

  # Runtime headers
  runtime/casadi_runtime.hpp
)

add_executable(casadi-cli casadi_cli.cpp)
target_link_libraries(casadi-cli casadi)
file(RELATIVE_PATH TREL_BIN_PREFIX "${CMAKE_INSTALL_PREFIX}" "${BIN_PREFIX}")
install(TARGETS casadi-cli
  RUNTIME DESTINATION ${TREL_BIN_PREFIX}
)

# Add the runtime sources to internal (installed separately)
set(RUNTIME_SRC)
foreach(FILE ${CASADI_RUNTIME_SRC})
  set(CASADI_INTERNAL ${CASADI_INTERNAL} runtime/${FILE})
  set(RUNTIME_SRC ${RUNTIME_SRC} runtime/${FILE})
endforeach()
set(RUNTIME_SRC ${RUNTIME_SRC} runtime/casadi_to_mex.hpp runtime/casadi_from_mex.hpp runtime/casadi_fmu.hpp)

string (REPLACE ";" "$<SEMICOLON>" ESCAPED_RUNTIME_SRC "${RUNTIME_SRC}")

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/runtime/casadi_runtime_str.h"
  COMMAND  ${CMAKE_COMMAND} -D OUTPUT="${CMAKE_CURRENT_BINARY_DIR}/runtime/casadi_runtime_str.h" -D SOURCES="${ESCAPED_RUNTIME_SRC}" -P "${CMAKE_CURRENT_SOURCE_DIR}/../generate_runtime.cmake"
  DEPENDS ${RUNTIME_SRC}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

set(CASADI_INTERNAL ${CASADI_INTERNAL} "${CMAKE_CURRENT_BINARY_DIR}/runtime/casadi_runtime_str.h")

# Build static and/or shared
casadi_library(casadi ${CASADI_PUBLIC} ${CASADI_INTERNAL})

file(RELATIVE_PATH TREL_INCLUDE_PREFIX "${CMAKE_INSTALL_PREFIX}" "${INCLUDE_PREFIX}")

target_include_directories(casadi
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/runtime>
    $<INSTALL_INTERFACE:${TREL_INCLUDE_PREFIX}>
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../> # <casadi.casadi.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../../> # <casadi/core/casadi_export.h>
)

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

if(CMAKE_CROSS_COMPILING AND MINGW)
  add_custom_command(TARGET casadi POST_BUILD
    COMMAND $ENV{DLLTOOL} --dllname $<TARGET_FILE:casadi> --def "${CMAKE_CURRENT_SOURCE_DIR}/casadi.def" --output-lib casadi.lib
  )
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/casadi.lib"
    DESTINATION ${DEF_LIB_PREFIX}
  )
endif()

install(FILES ${CASADI_PUBLIC} ${CMAKE_CURRENT_BINARY_DIR}/casadi_export.h
  DESTINATION ${INCLUDE_PREFIX}/casadi/core
)
if(INSTALL_INTERNAL_HEADERS)
  install(FILES ${CASADI_INTERNAL}
    DESTINATION ${INCLUDE_PREFIX}/casadi/core
  )
endif()

target_compile_definitions(casadi PUBLIC CASADI_SNPRINTF=snprintf)

if(WITH_THREADSAFE_SYMBOLICS)
  target_compile_definitions(casadi PUBLIC CASADI_WITH_THREADSAFE_SYMBOLICS)
endif()

if(WITH_THREAD)
  target_compile_definitions(casadi PUBLIC CASADI_WITH_THREAD)
endif()

if(MSVC)
  target_compile_options(casadi PRIVATE /bigobj)
endif()

get_target_property( compiler_definitions casadi COMPILE_DEFINITIONS )
foreach(DEF ${compiler_definitions})
  list(APPEND EXTRA_CXX_FLAGS_FROM_DEFS_LIST "-D${DEF}")
endforeach()
set(EXTRA_CXX_FLAGS_FROM_DEFS_LIST ${EXTRA_CXX_FLAGS_FROM_DEFS_LIST} PARENT_SCOPE)