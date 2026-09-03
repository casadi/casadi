# PATCH_COMMAND for blasfeo-external; runs with the working directory set to <SOURCE_DIR>.
#
# 1. blasfeo_alloc_guard.patch: the GENERIC lib4 kernels load full 4-wide blocks and mask
#    only on the write side, so they read up to one block past a dvec/dmat tail.
#
# 2. Rename the library to lib${CASADI_TP_SHARED_PREFIX}blasfeo. Upstream's SONAME is
#    libblasfeo.so.0, which every other blasfeo build also carries, so ld.so is free to
#    satisfy our libfatrop from a foreign copy loaded earlier in the process (acados ships
#    one). That is fatal because we build with EXT_DEP=OFF, which shrinks the public
#    blasfeo_timer typedef from a 32-byte struct to an int: a foreign blasfeo_tic/toc then
#    writes 32 bytes into fatrop's 16-byte Timer. See casadi/casadi#4395.
#
# Both steps are idempotent, so a re-run over an already-patched source tree is harmless.

if(NOT CASADI_TP_SHARED_PREFIX)
  message(FATAL_ERROR "patch_blasfeo.cmake: CASADI_TP_SHARED_PREFIX is not set")
endif()
if(NOT ALLOC_GUARD_PATCH)
  message(FATAL_ERROR "patch_blasfeo.cmake: ALLOC_GUARD_PATCH is not set")
endif()

execute_process(
  COMMAND git apply --ignore-whitespace "${ALLOC_GUARD_PATCH}"
  RESULT_VARIABLE PATCH_RESULT
  ERROR_VARIABLE PATCH_ERROR)
if(NOT PATCH_RESULT EQUAL 0)
  message(STATUS "blasfeo_alloc_guard.patch already applied or failed: ${PATCH_ERROR}")
endif()

file(READ CMakeLists.txt INPUT_TEXT)
if(NOT INPUT_TEXT MATCHES "CASADI_TP_RENAME")
  file(APPEND CMakeLists.txt
    "\n# CASADI_TP_RENAME: keep CasADi's private blasfeo off the libblasfeo SONAME.\n"
    "set_target_properties(blasfeo PROPERTIES OUTPUT_NAME \"${CASADI_TP_SHARED_PREFIX}blasfeo\")\n")
endif()
