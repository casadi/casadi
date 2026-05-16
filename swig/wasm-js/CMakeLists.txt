cmake_minimum_required(VERSION 3.10.2)

if(NOT EMSCRIPTEN)
  message(WARNING "WITH_WASM_JS=ON but build is not Emscripten; skipping. "
                  "Use 'emcmake cmake ..' to configure for WebAssembly.")
  return()
endif()

# Patched-SWIG -wasm-js capability check.
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

set(_iface     ${PROJECT_SOURCE_DIR}/swig/casadi.i)
set(_swig_cpp  ${CMAKE_CURRENT_BINARY_DIR}/casadi.cpp)
set(_swig_js   ${CMAKE_CURRENT_BINARY_DIR}/casadi.js)
set(_swig_dts  ${CMAKE_CURRENT_BINARY_DIR}/casadi.d.ts)
set(_swig_exp  ${CMAKE_CURRENT_BINARY_DIR}/casadi.cpp.exports)

add_custom_command(
  OUTPUT  ${_swig_cpp} ${_swig_js} ${_swig_dts} ${_swig_exp}
  COMMAND ${SWIG_EXECUTABLE} -wasm-js -stubs -c++
          -I${PROJECT_SOURCE_DIR}
          -I${PROJECT_SOURCE_DIR}/swig
          -o ${_swig_cpp} ${_iface}
  DEPENDS ${_iface}
  COMMENT "swig -wasm-js -stubs -> casadi.cpp + casadi.js + casadi.d.ts + exports"
  VERBATIM
)

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

# Build up the --preload-file list so SIDE_MODULE plugins are mounted in
# MEMFS at the path where casadi_os.cpp::open_shared_library expects them.
# `CASADI_SHARED_LIBRARY_PREFIX/SUFFIX` are baked into casadi_meta as
# "lib"/".so" regardless of host.  Mount each plugin at the root.
set(_preload_files "")

# ----------------------------------------------------------------------------
# Per-plugin SIDE_MODULE recipe
# ----------------------------------------------------------------------------
# wasmjs_side_plugin(<infix> <pname>) builds
# `lib${CMAKE_SHARED_LIBRARY_PREFIX}casadi_<infix>_<pname>.so` from the
# already-built static plugin archive (`casadi_<infix>_<pname>`),
# whole-archive'd so the `casadi_register_<infix>_<pname>` symbol
# survives DCE.  Each gets one extra .a (passed in as ${ARGN}) appended
# for solver-side third-party dependencies (e.g. libhighs.a).
# ----------------------------------------------------------------------------

function(wasmjs_side_plugin infix pname)
  set(_plugin_lib casadi_${infix}_${pname})
  set(_so_name    libcasadi_${infix}_${pname}.so)
  set(_so_path    ${CMAKE_CURRENT_BINARY_DIR}/${_so_name})
  # The casadi_plugin() macro emits its archive into ${CMAKE_BINARY_DIR}/lib.
  set(_plugin_a   ${CMAKE_BINARY_DIR}/lib/lib${_plugin_lib}.a)

  # Stub source so add_executable has something to compile.  Real content
  # is the --whole-archive of the static plugin lib at link time.
  set(_stub ${CMAKE_CURRENT_BINARY_DIR}/${_plugin_lib}_side_stub.c)
  file(WRITE ${_stub} "// Auto-generated SIDE_MODULE stub for ${_plugin_lib}.\n")

  add_executable(${_plugin_lib}_side ${_stub})
  add_dependencies(${_plugin_lib}_side ${_plugin_lib})
  target_compile_options(${_plugin_lib}_side PRIVATE -fPIC -fwasm-exceptions)

  set(_link_extra "")
  foreach(_extra ${ARGN})
    set(_link_extra "${_link_extra} ${_extra}")
  endforeach()

  # NOTE: $<TARGET_FILE:...> is NOT expanded inside the LINK_FLAGS
  # string property -- CMake only resolves generator expressions in a
  # handful of dedicated properties / commands (target_link_libraries,
  # target_link_options, etc.).  Hardcode the archive path that
  # casadi_plugin() emits to.  The add_dependencies() above still
  # forces the static lib to build first.
  # -Wl,--export-all: SIDE_MODULE=1 should auto-export everything but
  # in practice (emsdk circa 3.1.x+) DCE drops unreferenced exports
  # including each plugin's `casadi_register_<infix>_<pname>` symbol,
  # which is the entry point Plugin::load_plugin dlopens for.  Forcing
  # all symbols out preserves the registration symbols plus anything
  # the main module needs at dlopen time.
  set_target_properties(${_plugin_lib}_side PROPERTIES
    OUTPUT_NAME "casadi_${infix}_${pname}"
    PREFIX "lib"
    SUFFIX ".so"
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    LINK_FLAGS
      "-fPIC -fwasm-exceptions -s SIDE_MODULE=1 \
-Wl,--whole-archive ${_plugin_a} -Wl,--no-whole-archive \
-Wl,--export-all \
${_link_extra}"
  )
  # Propagate the .so out via a parent-scope variable so the main module
  # can --preload-file it.
  set(_preload_files
      "${_preload_files} --preload-file ${_so_path}@/${_so_name}"
      PARENT_SCOPE)
  add_dependencies(casadi_wasm ${_plugin_lib}_side)
endfunction()

# ----------------------------------------------------------------------------
# Plugins to ship as SIDE_MODULE .so's
# ----------------------------------------------------------------------------

# Built-in plugins implemented entirely in CasADi's own runtime (no
# third-party deps).  Each is one line: their static archive is just
# the plugin .cpp + its meta TU, and at runtime the side .so resolves
# everything else against the MAIN_MODULE.  Adds a few hundred KB
# total to MEMFS but unlocks a much larger surface in tests.
wasmjs_side_plugin(integrator rk)
wasmjs_side_plugin(integrator collocation)
wasmjs_side_plugin(conic       qrqp)
wasmjs_side_plugin(conic       ipqp)
wasmjs_side_plugin(conic       nlpsol)        # qp_to_nlp wrapper
wasmjs_side_plugin(nlpsol      sqpmethod)
wasmjs_side_plugin(nlpsol      qrsqp)
wasmjs_side_plugin(nlpsol      scpgen)
wasmjs_side_plugin(nlpsol      feasiblesqpmethod)
wasmjs_side_plugin(linsol      qr)
wasmjs_side_plugin(linsol      ldl)
wasmjs_side_plugin(linsol      tridiag)
wasmjs_side_plugin(linsol      lsqr)
wasmjs_side_plugin(linsol      symbolicqr)
wasmjs_side_plugin(interpolant linear)
wasmjs_side_plugin(interpolant bspline)
wasmjs_side_plugin(rootfinder  newton)
wasmjs_side_plugin(rootfinder  fast_newton)
wasmjs_side_plugin(rootfinder  nlpsol)        # implicit_to_nlp wrapper

# HiGHS Conic plugin: requires the cross-compiled libhighs.a archive.
# The parent's `highs` IMPORTED target points at a non-existent .so for
# the wasm build (ExternalProject_Add actually drops the .a there); use
# the .a's absolute path directly.
if(WITH_HIGHS)
  wasmjs_side_plugin(conic highs
    ${CMAKE_BINARY_DIR}/external_projects/lib/libhighs.a)
endif()

# DAQP Conic plugin: pure-C QP solver, no Fortran/BLAS dependencies.
if(WITH_DAQP)
  wasmjs_side_plugin(conic daqp
    ${CMAKE_BINARY_DIR}/external_projects/lib/libdaqp.a
    ${CMAKE_BINARY_DIR}/external_projects/lib/libdaqpstat.a)
endif()

# qpOASES Conic plugin: needs flang-built BLAS+LAPACK (dgemm_, dpotrf_, ...)
# plus the in-tree libcasadi_qpoases.a (qpOASES sources, patched to pass
# Fortran-77 hidden character-length args on EMSCRIPTEN).  Also pull in the
# Fortran runtime for the flang-built archives.
if(WITH_QPOASES)
  wasmjs_side_plugin(conic qpoases
    ${CMAKE_BINARY_DIR}/lib/libcasadi_qpoases.a
    ${CMAKE_BINARY_DIR}/external_projects/lib/liblapack.a
    ${CMAKE_BINARY_DIR}/external_projects/lib/libblas.a
    /opt/flang/wasm/lib/libFortranRuntime.a)
endif()

# Sundials integrators (cvodes, idas) and rootfinder (kinsol): pure-C,
# self-contained in the side module archive.
if(WITH_SUNDIALS)
  wasmjs_side_plugin(integrator cvodes)
  wasmjs_side_plugin(integrator idas)
  wasmjs_side_plugin(rootfinder kinsol)
endif()

# CSparse linsol: pure-C, self-contained.
if(WITH_CSPARSE)
  wasmjs_side_plugin(linsol csparse)
endif()

# LAPACK linsol plugins (lapacklu, lapackqr): need flang-built BLAS+LAPACK
# plus Fortran runtime.  Headers have been patched with EMSCRIPTEN-only
# Fortran-77 hidden character-length args.
if(WITH_LAPACK)
  wasmjs_side_plugin(linsol lapacklu
    ${CMAKE_BINARY_DIR}/external_projects/lib/liblapack.a
    ${CMAKE_BINARY_DIR}/external_projects/lib/libblas.a
    /opt/flang/wasm/lib/libFortranRuntime.a)
  wasmjs_side_plugin(linsol lapackqr
    ${CMAKE_BINARY_DIR}/external_projects/lib/liblapack.a
    ${CMAKE_BINARY_DIR}/external_projects/lib/libblas.a
    /opt/flang/wasm/lib/libFortranRuntime.a)
endif()

# Blocksqp Nlpsol plugin: uses qpoases internally for the inner QP solve,
# so needs the same BLAS+Fortran chain as qpoases.  Set qpsol_options/
# linsol to a wasm-available linsol (e.g. ldl) when constructing.
if(WITH_BLOCKSQP)
  wasmjs_side_plugin(nlpsol blocksqp
    ${CMAKE_BINARY_DIR}/lib/libcasadi_qpoases.a
    ${CMAKE_BINARY_DIR}/lib/libcasadi_conic_qpoases.a
    ${CMAKE_BINARY_DIR}/external_projects/lib/liblapack.a
    ${CMAKE_BINARY_DIR}/external_projects/lib/libblas.a
    /opt/flang/wasm/lib/libFortranRuntime.a)
endif()

# HPIPM Conic plugin: registers OK but requires OCP-block-structured QP
# inputs; no trivial smoke test possible.  Built but not currently
# exercised by test/javascript/*.
if(WITH_HPIPM)
  wasmjs_side_plugin(conic hpipm
    ${CMAKE_BINARY_DIR}/external_projects/lib/libhpipm.a
    ${CMAKE_BINARY_DIR}/external_projects/lib/libblasfeo.a)
endif()

# SuperSCS Conic plugin: registers OK but the SuperSCS sources call BLAS
# without the Fortran-77 hidden character-length args, triggering wasm-ld
# signature-mismatch -> unreachable on every BLAS-using call site.  Built
# but not operational until SuperSCS sources are patched the same way the
# vendored qpOASES was.
if(WITH_SUPERSCS)
  wasmjs_side_plugin(conic superscs
    ${CMAKE_BINARY_DIR}/external_projects/lib/libsuperscs.a
    ${CMAKE_BINARY_DIR}/external_projects/lib/libindirect.a
    ${CMAKE_BINARY_DIR}/external_projects/lib/liblinsys.a
    ${CMAKE_BINARY_DIR}/external_projects/lib/liblapack.a
    ${CMAKE_BINARY_DIR}/external_projects/lib/libblas.a
    /opt/flang/wasm/lib/libFortranRuntime.a)
endif()

# IPOPT Nlpsol plugin: pull in the full Fortran cascade.  libipopt.a +
# libcoinmumps.a + libcoinmetis.a from the COIN-OR/MUMPS/METIS chain,
# plus liblapack.a + libblas.a from Netlib reference LAPACK (built via
# r-wasm/flang-wasm), plus libFortranRuntime.a + the omp_lib stub
# object.  Order matters for the static-link symbol resolution -- ipopt
# pulls mumps, mumps pulls metis + lapack, lapack/blas pull Fortran
# runtime + omp stub.
if(WITH_IPOPT)
  set(_ipopt_extra
    ${CMAKE_BINARY_DIR}/external_projects/lib/libipopt.a
    ${CMAKE_BINARY_DIR}/external_projects/lib/libcoinmumps.a
    ${CMAKE_BINARY_DIR}/external_projects/lib/libcoinmetis.a
    ${CMAKE_BINARY_DIR}/external_projects/lib/liblapack.a
    ${CMAKE_BINARY_DIR}/external_projects/lib/libblas.a
    # libFortranRuntime.a from r-wasm/flang-wasm with our local
    # DescriptorModel.h patch baked in (see build-wasm-flang/
    # upstream-bug/bug_report.md).  The casadi-flang-wasm Docker image
    # replaces the shipped /opt/flang/wasm/lib/libFortranRuntime.a with
    # the patched build at image-build time.
    /opt/flang/wasm/lib/libFortranRuntime.a
    /opt/mumps_common_stub/mumps_common_stub.o)
  wasmjs_side_plugin(nlpsol ipopt ${_ipopt_extra})
endif()

# ----------------------------------------------------------------------------
# Main module link flags
# ----------------------------------------------------------------------------
# - MAIN_MODULE=1 gives us the dynamic linker at runtime (real dlopen).
# - The preload-file list mounts every plugin .so into MEMFS so the
#   first `Integrator("rk", ...)` / `qpsol("highs", ...)` call finds the
#   shared library on disk via the standard plugin search path.

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
-s EXPORTED_RUNTIME_METHODS=['UTF8ToString','stringToUTF8','lengthBytesUTF8','HEAPF64'] \
-s EXPORTED_FUNCTIONS=@${_swig_exp}\
${_preload_files}"
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
