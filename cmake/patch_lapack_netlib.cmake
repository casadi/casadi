# Patch Netlib reference LAPACK for wasm-emscripten builds via r-wasm flang.
#
# 1. LAPACKE/include/CMakeLists.txt invokes FortranCInterface_VERIFY(),
#    which tries to LINK a Fortran+C test executable.  Under
#    wasm32-unknown-emscripten this fails because the link line is
#    missing libFortranRuntime.a (the r-wasm flang runtime) and we
#    can't reach into Netlib's cmake to add it.  Remove the VERIFY
#    call; FortranCInterface_HEADER still works and produces a
#    correct `lapacke_mangling.h` using CMake's knowledge of the
#    flang mangling convention (lower + trailing underscore).
file(READ LAPACKE/include/CMakeLists.txt _TXT)
string(REPLACE
  "FortranCInterface_VERIFY()"
  "# patched-out by casadi cmake/patch_lapack_netlib.cmake (wasm build)"
  _TXT "${_TXT}")
file(WRITE LAPACKE/include/CMakeLists.txt "${_TXT}")

# 2. CBLAS/CMakeLists.txt has the same FortranCInterface_VERIFY() call.
if(EXISTS CBLAS/CMakeLists.txt)
  file(READ CBLAS/CMakeLists.txt _TXT)
  string(REPLACE
    "FortranCInterface_VERIFY()"
    "# patched-out by casadi cmake/patch_lapack_netlib.cmake (wasm build)"
    _TXT "${_TXT}")
  file(WRITE CBLAS/CMakeLists.txt "${_TXT}")
endif()
