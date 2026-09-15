cmake_minimum_required(VERSION 3.10.2)

# Three-stage build, mirroring swig/python and swig/matlab:
#   SWIG_EXPORT  -- run SWIG only; emit the wrapper (casadi.cpp/.js/.d.ts/
#                   .exports) into target/source for a downstream consumer.
#                   Runs in the SWIG codegen image (no Emscripten needed).
#   SWIG_IMPORT  -- compile the pre-generated wrapper against an already-built
#                   core; SWIG is NOT required (the wasm toolchain image ships
#                   no SWIG).  Emscripten.
#   neither      -- local all-in-one: run SWIG inline, then build.  Emscripten.

# Emscripten is required to *build* (IMPORT / local), but not to merely
# *generate* the wrapper sources (EXPORT).
if(NOT EMSCRIPTEN AND NOT SWIG_EXPORT)
  message(WARNING "WITH_WASM_JS=ON but build is not Emscripten; skipping. "
                  "Use 'emcmake cmake ..' to configure for WebAssembly.")
  return()
endif()

set(_iface     ${PROJECT_SOURCE_DIR}/swig/casadi.i)
# Pre-generated-wrapper surface a SWIG_IMPORT build consumes (committed/
# shipped by the SWIG_EXPORT producer), same layout as swig/python/target.
set(_wrap_dir  ${CMAKE_CURRENT_SOURCE_DIR}/target/source)

if(SWIG_IMPORT)
  message(STATUS "SWIG_IMPORT: compiling pre-generated wasm-js wrapper; SWIG not invoked")
  set(_swig_cpp  ${_wrap_dir}/casadi.cpp)
  set(_swig_js   ${_wrap_dir}/casadi.js)
  set(_swig_dts  ${_wrap_dir}/casadi.d.ts)
  set(_swig_exp  ${_wrap_dir}/casadi.cpp.exports)
  set_source_files_properties(${_swig_cpp} PROPERTIES GENERATED ON)
else()
  # Patched-SWIG -wasm-js capability check (only when we actually run SWIG).
  execute_process(
    COMMAND ${SWIG_EXECUTABLE} -help
    OUTPUT_VARIABLE _swig_help
    ERROR_VARIABLE  _swig_help_err
  )
  string(FIND "${_swig_help}${_swig_help_err}" "-wasm-js" _has_wasmjs)
  if(_has_wasmjs EQUAL -1)
    message(FATAL_ERROR
      "SWIG at ${SWIG_EXECUTABLE} lacks -wasm-js. Use the casadi-patched "
      "SWIG fork (see ~/programs/swig).")
  endif()

  if(SWIG_EXPORT)
    set(_out_dir ${_wrap_dir})
    file(MAKE_DIRECTORY ${_out_dir})
  else()
    set(_out_dir ${CMAKE_CURRENT_BINARY_DIR})
  endif()
  set(_swig_cpp  ${_out_dir}/casadi.cpp)
  set(_swig_js   ${_out_dir}/casadi.js)
  set(_swig_dts  ${_out_dir}/casadi.d.ts)
  set(_swig_exp  ${_out_dir}/casadi.cpp.exports)

  add_custom_command(
    OUTPUT  ${_swig_cpp} ${_swig_js} ${_swig_dts} ${_swig_exp}
    COMMAND ${SWIG_EXECUTABLE} -wasm-js -stubs -c++
            -I${PROJECT_SOURCE_DIR}
            -I${PROJECT_SOURCE_DIR}/swig
            -o ${_swig_cpp} ${_iface}
    # casadi.i's %insert("js") wraps load_<type>() to be async (fetch the
    # plugin .so into MEMFS, then dlopen).  The C++ signature is void, so
    # SWIG stubs them `: void`; retype to `: Promise<void>` to match.
    COMMAND ${CMAKE_COMMAND} -E env sed -i -E
            "s/(export function load_(nlpsol|conic|linsol|integrator|rootfinder|interpolant|expm|dple|blas|filesystem)\\(name: string\\)): void;/\\1: Promise<void>;/"
            ${_swig_dts}
    COMMAND ${CMAKE_COMMAND} -E env sed -i -E
            "s/(load_plugin\\(name: string\\)): void;/\\1: Promise<void>;/"
            ${_swig_dts}
    DEPENDS ${_iface} ${CMAKE_CURRENT_SOURCE_DIR}/load_plugins.js
    COMMENT "swig -wasm-js -stubs -> casadi.cpp + casadi.js + casadi.d.ts + exports"
    VERBATIM
  )

  if(SWIG_EXPORT)
    # Codegen-only target (the SWIG job builds this; no Emscripten core).
    add_custom_target(wasm_source DEPENDS ${_swig_cpp} ${_swig_js} ${_swig_dts} ${_swig_exp})
  endif()
endif()

# Under SWIG_EXPORT we only emit sources -- the main module + side-module
# plugins need the (Emscripten) core, which the codegen image doesn't build.
if(NOT SWIG_EXPORT)

# ----------------------------------------------------------------------------
# casadi_wasm MAIN_MODULE
# ----------------------------------------------------------------------------
# Build as Emscripten MAIN_MODULE so casadi_os.cpp's dlopen path resolves
# SIDE_MODULE plugin .so's at runtime.  This preserves CasADi's plugin-
# first / dlopen'd shared-library model in WASM -- GPL Ipopt+MUMPS,
# commercial CPLEX/Gurobi plugins, etc. are NEVER statically embedded in
# libcasadi (licensing nonstarter).
# ----------------------------------------------------------------------------

add_executable(casadi_wasm ${_swig_cpp})
add_dependencies(casadi_wasm casadi)
target_include_directories(casadi_wasm PRIVATE
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_BINARY_DIR}
)
target_link_libraries(casadi_wasm PRIVATE casadi)
target_compile_options(casadi_wasm PRIVATE -fwasm-exceptions -fPIC)

# Record the side modules that the asynchronous loader will fetch into MEMFS.
include(${CMAKE_CURRENT_SOURCE_DIR}/plugins.cmake)

# Link real targets so CMake carries dependencies and tracks archive changes.
function(wasmjs_side_plugin infix pname)
  set(_plugin casadi_${infix}_${pname})
  if(NOT TARGET ${_plugin})
    message(FATAL_ERROR "Missing configured plugin target ${_plugin}")
  endif()
  set(_side ${_plugin}_side)
  set(_stub ${wasm_dir}/${_side}_stub.cpp)
  file(WRITE ${_stub} "// Side module entry points come from the plugin archive.\n")
  add_executable(${_side} ${_stub})
  target_compile_options(${_side} PRIVATE -fPIC -fwasm-exceptions)
  target_link_libraries(${_side} PRIVATE
    "-Wl,--whole-archive" ${_plugin} "-Wl,--no-whole-archive")
  set_target_properties(${_side} PROPERTIES
    CASADI_WASM_SIDE_MODULE TRUE
    OUTPUT_NAME ${_plugin}
    PREFIX "lib"
    SUFFIX ".so"
    RUNTIME_OUTPUT_DIRECTORY ${wasm_dir}
    LINK_FLAGS "-fPIC -fwasm-exceptions -s SIDE_MODULE=1 -Wl,--export-all")
  wasmjs_register_plugin(${infix} ${pname})
  add_dependencies(casadi_wasm ${_side})
endfunction()

# Called from the root after casadi/ has registered all plugin targets.
function(wasmjs_finalize_plugins wasm_dir)
  # flang needs explicit storage for MUMPS's libseq COMMON block.
  if(WITH_IPOPT OR WITH_MUMPS)
    set(_mumps_stub_obj ${wasm_dir}/mumps_common_stub.o)
    add_custom_command(OUTPUT ${_mumps_stub_obj}
      COMMAND ${CMAKE_Fortran_COMPILER} -fPIC -c
        ${PROJECT_SOURCE_DIR}/swig/wasm-js/mumps_common_stub.f90 -o ${_mumps_stub_obj}
      DEPENDS ${PROJECT_SOURCE_DIR}/swig/wasm-js/mumps_common_stub.f90
      VERBATIM)
    add_custom_target(mumps_common_stub DEPENDS ${_mumps_stub_obj})
    set_property(TARGET mumps APPEND PROPERTY INTERFACE_LINK_LIBRARIES ${_mumps_stub_obj})
    add_dependencies(mumps mumps_common_stub)
  endif()

  # Built-in plugins have no additional platform requirements.
  foreach(_entry
      "integrator::rk" "integrator::collocation"
      "conic::qrqp" "conic::ipqp" "conic::nlpsol"
      "nlpsol::sqpmethod" "nlpsol::qrsqp" "nlpsol::scpgen" "nlpsol::feasiblesqpmethod"
      "linsol::qr" "linsol::ldl" "linsol::tridiag" "linsol::lsqr" "linsol::symbolicqr"
      "interpolant::linear" "interpolant::bspline"
      "rootfinder::newton" "rootfinder::fast_newton" "rootfinder::nlpsol"
      "rootfinder::bisection")
    string(REPLACE "::" ";" _parts "${_entry}")
    wasmjs_side_plugin(${_parts})
  endforeach()

  # Known external plugins use the same link recipe, including transitive dependencies.
  foreach(_entry
      "conic::highs" "conic::daqp" "conic::qpoases" "conic::hpipm" "conic::osqp"
      "conic::fatrop" "nlpsol::fatrop" "nlpsol::ipopt" "nlpsol::blocksqp"
      "integrator::cvodes" "integrator::idas" "rootfinder::kinsol"
      "linsol::csparse" "linsol::csparsecholesky" "linsol::lapacklu" "linsol::lapackqr"
      "linsol::mumps" "blas::classic" "blas::blasfeo" "xmlfile::tinyxml")
    string(REPLACE "::" "_" _target "${_entry}")
    if(TARGET casadi_${_target})
      string(REPLACE "::" ";" _parts "${_entry}")
      wasmjs_side_plugin(${_parts})
    endif()
  endforeach()

  # Every enabled plugin must ship or have a documented browser exclusion.
  wasmjs_exclude_plugin(importer shell "Requires spawning a native C compiler")
  wasmjs_exclude_plugin(conic superscs "Requires fixes for the Fortran CHARACTER-length ABI")
  wasmjs_write_manifest(${wasm_dir}/plugins.json)

endfunction()

# ----------------------------------------------------------------------------
# Main module link flags
# ----------------------------------------------------------------------------
# - MAIN_MODULE=1 gives us the dynamic linker at runtime (real dlopen).
# - Plugins are fetched into MEMFS on demand by the JavaScript loader.

set_target_properties(casadi_wasm PROPERTIES
  OUTPUT_NAME "casadi_wasm"
  SUFFIX ".js"
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  LINK_FLAGS "-fwasm-exceptions \
-lembind \
-fPIC -s MAIN_MODULE=1 \
-Wl,--export-all \
-s MODULARIZE=1 -s EXPORT_NAME=createWasm \
-s ENVIRONMENT=web,worker,node \
-s ALLOW_MEMORY_GROWTH=1 -s ASSERTIONS=1 \
-s FORCE_FILESYSTEM=1 \
-s EXPORTED_RUNTIME_METHODS=['UTF8ToString','stringToUTF8','lengthBytesUTF8','HEAPF64','FS'] \
-s EXPORTED_FUNCTIONS=@${_swig_exp}"
)
# -Wl,--export-all on the MAIN_MODULE: same rationale as the side
# modules above.  MAIN_MODULE=1 is supposed to export everything for
# side-module dlopen resolution, but DCE silently drops static-data
# symbols (e.g. `casadi::Conic::options_`) when nothing in the wasm
# wrapper code references them.  Plugin dlopens then fail with
# `Aborted(Assertion failed: undefined symbol '_ZN6casadi5Conic8options_E'.
# perhaps a side module was not linked in?...)`.  Forcing all symbols
# out exports the whole libcasadi.a surface.

# Stage demo assets next to the linked module so a static server can serve
# them directly.  Runtime + typecheck tests live in test/javascript/ and
# load the auto-generated casadi.js (no hand-rolled JS layer required:
# wasm_js.cxx emits class-return wrapping so methods return proxy
# instances instead of raw void* pointers).
foreach(_asset index.html)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/${_asset}
    ${CMAKE_CURRENT_BINARY_DIR}/${_asset}
    COPYONLY)
endforeach()

endif() # NOT SWIG_EXPORT
