cmake_minimum_required(VERSION 2.8.6)

    if(WIN32)
    set(BLASFEO_FLAGS -DOS_WINDOWS)
    elseif(APPLE)
     set(BLASFEO_FLAGS -DOS_MAC)
    else()
     set(BLASFEO_FLAGS -DOS_LINUX)
    endif()
    if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
casadi_external_library(casadi_blasfeo_X64_INTEL_HASWELL src/auxiliary/d_aux_lib4.c
src/auxiliary/d_aux_extern_depend_lib4.c
src/auxiliary/s_aux_lib4.c
src/auxiliary/s_aux_extern_depend_lib4.c
src/auxiliary/i_aux_extern_depend_lib4.c
src/auxiliary/avx2/kernel_dgetr_lib4.c
src/auxiliary/avx/kernel_dgecp_lib4.c
src/kernel/avx2/kernel_dgemm_4x4_lib4.S
src/kernel/avx2/kernel_dgemm_8x4_lib4.S
src/kernel/avx2/kernel_dgemm_12x4_lib4.S
src/kernel/avx2/kernel_dgemv_8_lib4.S
src/kernel/avx2/kernel_dgetrf_pivot_4_lib4.c
src/kernel/avx2/kernel_dsymv_6_lib4.S
src/kernel/avx/kernel_dgemv_4_lib4.S
src/kernel/c99/kernel_dsymv_4_lib4.c
src/kernel/c99/kernel_sgemm_4x4_lib4.c
src/kernel/c99/kernel_sgemv_4_lib4.c
src/blas/d_blas3_lib4.c
src/blas/d_lapack_lib4.c
src/blas/d_blas2_lib4.c
src/blas/d_blas1_lib4.c
src/blas/s_blas3_lib4.c
src/blas/s_lapack_lib4.c
src/blas/s_blas2_lib4.c)
set_target_properties(casadi_blasfeo_X64_INTEL_HASWELL PROPERTIES COMPILE_FLAGS "-O2 -fPIC -m64 -mavx2 -mfma -DTARGET_X64_INTEL_HASWELL -DLA_BLASFEO -DREF_BLAS_OPENBLAS ${BLASFEO_FLAGS}")

set_property(SOURCE src/kernel/avx2/kernel_dgemm_4x4_lib4.S PROPERTY LANGUAGE C)
set_property(SOURCE src/kernel/avx2/kernel_dgemm_8x4_lib4.S PROPERTY LANGUAGE C)
set_property(SOURCE src/kernel/avx2/kernel_dgemm_12x4_lib4.S PROPERTY LANGUAGE C)
set_property(SOURCE src/kernel/avx2/kernel_dgemv_8_lib4.S PROPERTY LANGUAGE C)
set_property(SOURCE src/kernel/avx2/kernel_dsymv_6_lib4.S PROPERTY LANGUAGE C)
set_property(SOURCE src/kernel/avx/kernel_dgemv_4_lib4.S PROPERTY LANGUAGE C)
endif()
if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
casadi_external_library(casadi_blasfeo_X64_INTEL_SANDY_BRIDGE src/auxiliary/d_aux_lib4.c
src/auxiliary/d_aux_extern_depend_lib4.c
src/auxiliary/s_aux_lib4.c
src/auxiliary/s_aux_extern_depend_lib4.c
src/auxiliary/i_aux_extern_depend_lib4.c
src/auxiliary/avx/kernel_dgecp_lib4.c
src/auxiliary/avx/kernel_dgetr_lib4.c
src/kernel/avx/kernel_dgemm_4x4_lib4.S
src/kernel/avx/kernel_dgemm_8x4_lib4.S
src/kernel/avx/kernel_dgemv_12_lib4.S
src/kernel/avx/kernel_dgemv_8_lib4.S
src/kernel/avx/kernel_dgemv_4_lib4.S
src/kernel/avx/kernel_dgetrf_pivot_4_lib4.c
src/kernel/avx/kernel_dsymv_6_lib4.S
src/kernel/c99/kernel_dsymv_4_lib4.c
src/kernel/c99/kernel_sgemm_4x4_lib4.c
src/kernel/c99/kernel_sgemv_4_lib4.c
src/blas/d_blas3_lib4.c
src/blas/d_lapack_lib4.c
src/blas/d_blas2_lib4.c
src/blas/d_blas1_lib4.c
src/blas/s_blas3_lib4.c
src/blas/s_lapack_lib4.c
src/blas/s_blas2_lib4.c)
set_target_properties(casadi_blasfeo_X64_INTEL_SANDY_BRIDGE PROPERTIES COMPILE_FLAGS "-O2 -fPIC -m64 -mavx -DTARGET_X64_INTEL_SANDY_BRIDGE -DLA_BLASFEO -DREF_BLAS_OPENBLAS ${BLASFEO_FLAGS}")

set_property(SOURCE src/kernel/avx/kernel_dgemm_4x4_lib4.S PROPERTY LANGUAGE C)
set_property(SOURCE src/kernel/avx/kernel_dgemm_8x4_lib4.S PROPERTY LANGUAGE C)
set_property(SOURCE src/kernel/avx/kernel_dgemv_12_lib4.S PROPERTY LANGUAGE C)
set_property(SOURCE src/kernel/avx/kernel_dgemv_8_lib4.S PROPERTY LANGUAGE C)
set_property(SOURCE src/kernel/avx/kernel_dgemv_4_lib4.S PROPERTY LANGUAGE C)
set_property(SOURCE src/kernel/avx/kernel_dsymv_6_lib4.S PROPERTY LANGUAGE C)
endif()
if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
casadi_external_library(casadi_blasfeo_X64_INTEL_CORE src/auxiliary/d_aux_lib4.c
src/auxiliary/d_aux_extern_depend_lib4.c
src/auxiliary/s_aux_lib4.c
src/auxiliary/s_aux_extern_depend_lib4.c
src/auxiliary/i_aux_extern_depend_lib4.c
src/auxiliary/c99/kernel_dgecp_lib4.c
src/auxiliary/c99/kernel_dgetr_lib4.c
src/kernel/sse3/kernel_dgemm_4x4_lib4.S
src/kernel/c99/kernel_dgemm_4x4_lib4.c
src/kernel/c99/kernel_dgemv_4_lib4.c
src/kernel/c99/kernel_dgetrf_pivot_4_lib4.c
src/kernel/c99/kernel_dsymv_4_lib4.c
src/kernel/c99/kernel_sgemm_4x4_lib4.c
src/kernel/c99/kernel_sgemv_4_lib4.c
src/blas/d_blas3_lib4.c
src/blas/d_lapack_lib4.c
src/blas/d_blas2_lib4.c
src/blas/d_blas1_lib4.c
src/blas/s_blas3_lib4.c
src/blas/s_lapack_lib4.c
src/blas/s_blas2_lib4.c)
set_target_properties(casadi_blasfeo_X64_INTEL_CORE PROPERTIES COMPILE_FLAGS "-O2 -fPIC -m64 -msse3 -DTARGET_X64_INTEL_CORE -DLA_BLASFEO -DREF_BLAS_OPENBLAS ${BLASFEO_FLAGS}")

set_property(SOURCE src/kernel/sse3/kernel_dgemm_4x4_lib4.S PROPERTY LANGUAGE C)
endif()
if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
casadi_external_library(casadi_blasfeo_X64_AMD_BULLDOZER src/auxiliary/d_aux_lib4.c
src/auxiliary/d_aux_extern_depend_lib4.c
src/auxiliary/s_aux_lib4.c
src/auxiliary/s_aux_extern_depend_lib4.c
src/auxiliary/i_aux_extern_depend_lib4.c
src/auxiliary/c99/kernel_dgecp_lib4.c
src/auxiliary/c99/kernel_dgetr_lib4.c
src/kernel/fma/kernel_dgemm_4x4_lib4.S
src/kernel/c99/kernel_dgemm_4x4_lib4.c
src/kernel/c99/kernel_dgemv_4_lib4.c
src/kernel/c99/kernel_dgetrf_pivot_4_lib4.c
src/kernel/c99/kernel_dsymv_4_lib4.c
src/kernel/c99/kernel_sgemm_4x4_lib4.c
src/kernel/c99/kernel_sgemv_4_lib4.c
src/blas/d_blas3_lib4.c
src/blas/d_lapack_lib4.c
src/blas/d_blas2_lib4.c
src/blas/d_blas1_lib4.c
src/blas/s_blas3_lib4.c
src/blas/s_lapack_lib4.c
src/blas/s_blas2_lib4.c)
set_target_properties(casadi_blasfeo_X64_AMD_BULLDOZER PROPERTIES COMPILE_FLAGS "-O2 -fPIC -m64 -mavx -mfma -DTARGET_X64_AMD_BULLDOZER -DLA_BLASFEO -DREF_BLAS_OPENBLAS ${BLASFEO_FLAGS}")

set_property(SOURCE src/kernel/fma/kernel_dgemm_4x4_lib4.S PROPERTY LANGUAGE C)
endif()
casadi_external_library(casadi_blasfeo_GENERIC src/auxiliary/d_aux_lib4.c
src/auxiliary/d_aux_extern_depend_lib4.c
src/auxiliary/s_aux_lib4.c
src/auxiliary/s_aux_extern_depend_lib4.c
src/auxiliary/i_aux_extern_depend_lib4.c
src/auxiliary/c99/kernel_dgecp_lib4.c
src/auxiliary/c99/kernel_dgetr_lib4.c
src/kernel/c99/kernel_dgemm_4x4_lib4.c
src/kernel/c99/kernel_dgemv_4_lib4.c
src/kernel/c99/kernel_dgetrf_pivot_4_lib4.c
src/kernel/c99/kernel_dsymv_4_lib4.c
src/kernel/c99/kernel_sgemm_4x4_lib4.c
src/kernel/c99/kernel_sgemv_4_lib4.c
src/blas/d_blas3_lib4.c
src/blas/d_lapack_lib4.c
src/blas/d_blas2_lib4.c
src/blas/d_blas1_lib4.c
src/blas/s_blas3_lib4.c
src/blas/s_lapack_lib4.c
src/blas/s_blas2_lib4.c)
set_target_properties(casadi_blasfeo_GENERIC PROPERTIES COMPILE_FLAGS "-O2 -fPIC -march=native -DTARGET_GENERIC -DLA_BLASFEO -DREF_BLAS_OPENBLAS ${BLASFEO_FLAGS}")

