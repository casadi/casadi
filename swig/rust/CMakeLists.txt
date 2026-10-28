cmake_minimum_required(VERSION 3.10.2)

# Three-stage build, as swig/julia:
#   SWIG_EXPORT  -- run SWIG only; emit casadiRUST_wrap.cxx + casadi.rs into gen/
#   SWIG_IMPORT  -- compile the pre-generated wrapper against an already-built core
#   neither      -- run SWIG inline, then build
# The crate (Cargo.toml here) builds against the directory holding casadi.rs and
# libcasadi_rust_wrap: CASADI_RUST_GEN=<build>/swig/rust.

set(_iface    ${PROJECT_SOURCE_DIR}/swig/casadi.i)
set(_wrap_dir ${CMAKE_CURRENT_SOURCE_DIR}/gen)

if(SWIG_IMPORT)
  message(STATUS "SWIG_IMPORT: compiling pre-generated Rust wrapper; SWIG not invoked")
  set(_out_dir ${_wrap_dir})
else()
  execute_process(COMMAND ${SWIG_EXECUTABLE} -help OUTPUT_VARIABLE _swig_help ERROR_VARIABLE _swig_help_err)
  string(FIND "${_swig_help}${_swig_help_err}" "-rust" _has_rust)
  if(_has_rust EQUAL -1)
    message(FATAL_ERROR "SWIG at ${SWIG_EXECUTABLE} lacks -rust. Use the casadi-patched SWIG (ghcr.io/casadi/ci-swig).")
  endif()
  if(SWIG_EXPORT)
    set(_out_dir ${_wrap_dir})
  else()
    set(_out_dir ${CMAKE_CURRENT_BINARY_DIR})
  endif()
  file(MAKE_DIRECTORY ${_out_dir})
endif()

set(_swig_cpp ${_out_dir}/casadiRUST_wrap.cxx)
set(_swig_rs  ${_out_dir}/casadi.rs)

if(NOT SWIG_IMPORT)
  configure_file(../swig_config.h.cmake "${CMAKE_CURRENT_BINARY_DIR}/swig_config.h" ESCAPE_QUOTES)
  add_custom_command(
    OUTPUT  ${_swig_cpp} ${_swig_rs}
    COMMAND ${SWIG_EXECUTABLE} -rust -c++ ${CMAKE_SWIG_FLAGS} -w524,503 ${MAYBE_WERROR}
            -I${CMAKE_CURRENT_BINARY_DIR}
            -I${PROJECT_SOURCE_DIR}
            -I${PROJECT_SOURCE_DIR}/swig
            -outdir ${_out_dir}
            -o ${_swig_cpp} ${_iface}
    DEPENDS ${_iface}
    COMMENT "swig -rust -> casadiRUST_wrap.cxx + casadi.rs"
    VERBATIM
  )
  if(SWIG_EXPORT)
    add_custom_target(rust_source DEPENDS ${_swig_cpp} ${_swig_rs})
    return()
  endif()
endif()

add_library(casadi_rust SHARED ${_swig_cpp})
set_target_properties(casadi_rust PROPERTIES
  PREFIX "lib" OUTPUT_NAME "casadi_rust_wrap"
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  CXX_STANDARD 17)
target_include_directories(casadi_rust PRIVATE ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR})
target_compile_definitions(casadi_rust PRIVATE WITH_DEPRECATED_FEATURES)
target_link_libraries(casadi_rust casadi)

# casadi.rs next to the library: the crate's CASADI_RUST_GEN
add_custom_command(TARGET casadi_rust POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${_swig_rs} ${CMAKE_CURRENT_BINARY_DIR}/casadi.rs
  COMMENT "casadi.rs -> ${CMAKE_CURRENT_BINARY_DIR}")

install(TARGETS casadi_rust DESTINATION ${LIB_PREFIX})
install(FILES ${_swig_rs} DESTINATION ${LIB_PREFIX}/rust)
