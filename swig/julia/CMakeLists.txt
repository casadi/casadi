cmake_minimum_required(VERSION 3.10.2)

# Three-stage build, mirroring swig/wasm-js:
#   SWIG_EXPORT  -- run SWIG only; emit casadi_wrap.cxx + casadi.jl into
#                   target/source (runs in the SWIG codegen image; no Julia
#                   needed -- generation does not touch julia.h).
#   SWIG_IMPORT  -- compile the pre-generated wrapper against an already-built
#                   core; SWIG not required.  Needs Julia headers.
#   neither      -- local all-in-one: run SWIG inline, then build.

set(_iface     ${PROJECT_SOURCE_DIR}/swig/casadi.i)
set(_wrap_dir  ${CMAKE_CURRENT_SOURCE_DIR}/target/source)

if(SWIG_IMPORT)
  message(STATUS "SWIG_IMPORT: compiling pre-generated julia wrapper; SWIG not invoked")
  set(_out_dir ${_wrap_dir})
else()
  find_package(SWIG REQUIRED)
  execute_process(COMMAND ${SWIG_EXECUTABLE} -help OUTPUT_VARIABLE _swig_help ERROR_VARIABLE _swig_help_err)
  string(FIND "${_swig_help}${_swig_help_err}" "-julia" _has_julia)
  if(_has_julia EQUAL -1)
    message(FATAL_ERROR
      "SWIG at ${SWIG_EXECUTABLE} lacks -julia. Use the casadi-patched "
      "SWIG (ghcr.io/casadi/ci-swig).")
  endif()
  if(SWIG_EXPORT)
    set(_out_dir ${_wrap_dir})
    file(MAKE_DIRECTORY ${_out_dir})
  else()
    set(_out_dir ${CMAKE_CURRENT_BINARY_DIR})
  endif()
endif()

set(_swig_cpp ${_out_dir}/casadi_wrap.cxx)
set(_swig_jl  ${_out_dir}/casadi.jl)

if(NOT SWIG_IMPORT)
  add_custom_command(
    OUTPUT  ${_swig_cpp} ${_swig_jl}
    COMMAND ${SWIG_EXECUTABLE} -julia -c++
            -I${PROJECT_SOURCE_DIR}
            -I${PROJECT_SOURCE_DIR}/swig
            -o ${_swig_cpp} ${_iface}
    DEPENDS ${_iface}
    COMMENT "swig -julia -> casadi_wrap.cxx + casadi.jl"
    VERBATIM
  )
  if(SWIG_EXPORT)
    # Codegen-only target (the SWIG job builds this; no Julia needed).
    add_custom_target(julia_source DEPENDS ${_swig_cpp} ${_swig_jl})
    return()
  endif()
endif()

# ---- compile stage: needs Julia headers (julia.h; libjulia is NOT linked,
# jl_* symbols resolve from the host process at dlopen) ----
if(NOT JULIA_INCLUDE_DIR)
  find_program(JULIA_EXECUTABLE julia)
  if(JULIA_EXECUTABLE)
    execute_process(
      COMMAND ${JULIA_EXECUTABLE} -e "print(abspath(Sys.BINDIR, \"..\", \"include\", \"julia\"))"
      OUTPUT_VARIABLE JULIA_INCLUDE_DIR)
  endif()
endif()
if(NOT JULIA_INCLUDE_DIR OR NOT EXISTS "${JULIA_INCLUDE_DIR}/julia.h")
  message(FATAL_ERROR "julia.h not found; set -DJULIA_INCLUDE_DIR=<julia>/include/julia")
endif()
message(STATUS "Julia headers: ${JULIA_INCLUDE_DIR}")

add_library(casadi_julia SHARED ${_swig_cpp})
set_target_properties(casadi_julia PROPERTIES
  PREFIX "lib" OUTPUT_NAME "casadi_wrap"
  CXX_STANDARD 17)
target_include_directories(casadi_julia PRIVATE
  ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} ${JULIA_INCLUDE_DIR})
target_compile_definitions(casadi_julia PRIVATE WITH_DEPRECATED_FEATURES)
target_link_libraries(casadi_julia casadi)

# Assemble the consumable module dir: casadi.jl + CasADi.jl + the wrapper .so
set(_jl_pkg ${CMAKE_CURRENT_BINARY_DIR}/CasADiNative)
add_custom_command(TARGET casadi_julia POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E make_directory ${_jl_pkg}
  COMMAND ${CMAKE_COMMAND} -E copy ${_swig_jl} ${_jl_pkg}/casadi.jl
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/CasADiNative.jl ${_jl_pkg}/CasADiNative.jl
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:casadi_julia> ${_jl_pkg}/libcasadi_wrap.so
  COMMENT "assemble Julia module dir ${_jl_pkg}")

install(TARGETS casadi_julia DESTINATION ${LIB_PREFIX})
install(FILES ${_swig_jl} ${CMAKE_CURRENT_SOURCE_DIR}/CasADiNative.jl DESTINATION ${LIB_PREFIX}/julia)
