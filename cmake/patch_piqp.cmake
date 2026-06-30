# Allow PIQP to configure with CMake 3.20 (the toolchain CMake on the mingw/mxe
# Windows cross-build is 3.20.1, while PIQP pins cmake_minimum_required 3.21).
#
# 3.21 is not actually required for the casadi build: the only >3.20 feature in
# PIQP's CMake is the $<TARGET_RUNTIME_DLLS> generator expression (added in
# CMake 3.21) inside the WIN32-only fix_test_dll macro, which is invoked solely
# by the test targets. We build with BUILD_TESTS=OFF, so it is never reached;
# the guard below additionally keeps configuration working if tests are enabled
# on CMake 3.20 (DLLs simply are not auto-copied next to the test binaries).

# Lower every cmake_minimum_required(VERSION 3.21) in the tree to 3.20.
file(GLOB_RECURSE PIQP_CMAKE_LISTS "${SRC}/*/CMakeLists.txt")
list(APPEND PIQP_CMAKE_LISTS "${SRC}/CMakeLists.txt")
foreach(f ${PIQP_CMAKE_LISTS})
  file(READ "${f}" txt)
  string(REGEX REPLACE
    "cmake_minimum_required\\(VERSION 3\\.21\\)"
    "cmake_minimum_required(VERSION 3.20)"
    txt "${txt}")
  file(WRITE "${f}" "${txt}")
endforeach()

# Guard the lone $<TARGET_RUNTIME_DLLS> usage behind a CMake 3.21 check.
set(PIQP_TOP "${SRC}/CMakeLists.txt")
file(READ "${PIQP_TOP}" txt)
string(REPLACE
  "    if (WIN32)\n        add_custom_command("
  "    if (WIN32 AND CMAKE_VERSION VERSION_GREATER_EQUAL \"3.21\")\n        add_custom_command("
  txt "${txt}")
file(WRITE "${PIQP_TOP}" "${txt}")

# Zero-initialise the bound-index arrays in {sparse,dense} Data::resize().
# PIQP fills h_l_idx/h_u_idx/x_l_idx/x_u_idx only for FINITE bounds; the tail
# (for +/-inf bounds, i.e. essentially every LP) stays uninitialised and is
# read in KKTSystem::solve. Harmless on glibc, but the garbage intermittently
# derails the interior-point solve on the mingw Windows build (flaky LP
# failures). Verified with valgrind: this drops 6 -> 0 uninitialised reads.
foreach(dt "sparse/data.tpp" "dense/data.tpp")
  set(f "${SRC}/include/piqp/${dt}")
  file(READ "${f}" txt)
  string(REPLACE
    "    h_l_idx.resize(m);\n    h_u_idx.resize(m);\n    x_l_idx.resize(n);\n    x_u_idx.resize(n);"
    "    h_l_idx.resize(m);\n    h_u_idx.resize(m);\n    x_l_idx.resize(n);\n    x_u_idx.resize(n);\n    h_l_idx.setZero(); h_u_idx.setZero(); x_l_idx.setZero(); x_u_idx.setZero();"
    txt "${txt}")
  file(WRITE "${f}" "${txt}")
endforeach()
