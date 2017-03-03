cmake_minimum_required(VERSION 2.8.6)

    if(WIN32)
    set(HPMPC_FLAGS -DOS_WINDOWS)
    elseif(APPLE)
     set(HPMPC_FLAGS -DOS_MAC)
    else()
     set(HPMPC_FLAGS -DOS_LINUX)
    endif()
    if()
    set(CMAKE_SHARED_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS} -Wl,--unresolved-symbols=report-all)
    endif()
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../blasfeo/src/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/dummy_include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/dummy_include/foo/bar)
if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
casadi_external_library(casadi_hpmpc_X64_AVX_GENERIC src/auxiliary/d_aux_lib4.c
src/auxiliary/d_aux_extern_depend_lib4.c
src/auxiliary/i_aux.c
src/lqcp_solvers/d_back_ric_rec_libstr.c
src/lqcp_solvers/d_part_cond_libstr.c
src/lqcp_solvers/d_tree_back_ric_rec_libstr.c
src/mpc_solvers/d_ip2_res_hard_libstr.c
src/mpc_solvers/d_tree_ip2_res_hard_libstr.c
src/mpc_solvers/d_res_ip_res_hard_libstr.c
src/mpc_solvers/d_tree_res_ip_res_hard_libstr.c
src/mpc_solvers/avx/d_aux_ip_hard_libstr.c
src/interfaces/c/fortran_order_interface_libstr.c)
SET_TARGET_PROPERTIES(casadi_hpmpc_X64_AVX_GENERIC PROPERTIES COMPILE_FLAGS "-DLA_BLASFEO -DTARGET_GENERIC -DTARGET_X64_AVX -O2 -fPIC -DBLASFEO -m64 -mavx -DTARGET_X64_AVX ${HPMPC_FLAGS}")

target_link_libraries(casadi_hpmpc_X64_AVX_GENERIC casadi_blasfeo_GENERIC m)
endif()
if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
casadi_external_library(casadi_hpmpc_X64_AVX2_X64_INTEL_HASWELL src/auxiliary/d_aux_lib4.c
src/auxiliary/d_aux_extern_depend_lib4.c
src/auxiliary/i_aux.c
src/lqcp_solvers/d_back_ric_rec_libstr.c
src/lqcp_solvers/d_part_cond_libstr.c
src/lqcp_solvers/d_tree_back_ric_rec_libstr.c
src/mpc_solvers/d_ip2_res_hard_libstr.c
src/mpc_solvers/d_tree_ip2_res_hard_libstr.c
src/mpc_solvers/d_res_ip_res_hard_libstr.c
src/mpc_solvers/d_tree_res_ip_res_hard_libstr.c
src/mpc_solvers/avx/d_aux_ip_hard_libstr.c
src/interfaces/c/fortran_order_interface_libstr.c)
SET_TARGET_PROPERTIES(casadi_hpmpc_X64_AVX2_X64_INTEL_HASWELL PROPERTIES COMPILE_FLAGS "-DLA_BLASFEO -DTARGET_X64_INTEL_HASWELL -DTARGET_X64_AVX2 -O2 -fPIC -DBLASFEO -m64 -mavx2 -mfma -DTARGET_X64_AVX2 ${HPMPC_FLAGS}")

target_link_libraries(casadi_hpmpc_X64_AVX2_X64_INTEL_HASWELL casadi_blasfeo_X64_INTEL_HASWELL m)
endif()
if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
casadi_external_library(casadi_hpmpc_X64_AVX_X64_INTEL_SANDY_BRIDGE src/auxiliary/d_aux_lib4.c
src/auxiliary/d_aux_extern_depend_lib4.c
src/auxiliary/i_aux.c
src/lqcp_solvers/d_back_ric_rec_libstr.c
src/lqcp_solvers/d_part_cond_libstr.c
src/lqcp_solvers/d_tree_back_ric_rec_libstr.c
src/mpc_solvers/d_ip2_res_hard_libstr.c
src/mpc_solvers/d_tree_ip2_res_hard_libstr.c
src/mpc_solvers/d_res_ip_res_hard_libstr.c
src/mpc_solvers/d_tree_res_ip_res_hard_libstr.c
src/mpc_solvers/avx/d_aux_ip_hard_libstr.c
src/interfaces/c/fortran_order_interface_libstr.c)
SET_TARGET_PROPERTIES(casadi_hpmpc_X64_AVX_X64_INTEL_SANDY_BRIDGE PROPERTIES COMPILE_FLAGS "-DLA_BLASFEO -DTARGET_X64_INTEL_SANDY_BRIDGE -DTARGET_X64_AVX -O2 -fPIC -DBLASFEO -m64 -mavx -DTARGET_X64_AVX ${HPMPC_FLAGS}")

target_link_libraries(casadi_hpmpc_X64_AVX_X64_INTEL_SANDY_BRIDGE casadi_blasfeo_X64_INTEL_SANDY_BRIDGE m)
endif()
if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
casadi_external_library(casadi_hpmpc_X64_SSE3_X64_INTEL_CORE src/auxiliary/d_aux_lib4.c
src/auxiliary/d_aux_extern_depend_lib4.c
src/auxiliary/i_aux.c
src/lqcp_solvers/d_back_ric_rec_libstr.c
src/lqcp_solvers/d_part_cond_libstr.c
src/lqcp_solvers/d_tree_back_ric_rec_libstr.c
src/mpc_solvers/d_ip2_res_hard_libstr.c
src/mpc_solvers/d_tree_ip2_res_hard_libstr.c
src/mpc_solvers/d_res_ip_res_hard_libstr.c
src/mpc_solvers/d_tree_res_ip_res_hard_libstr.c
src/mpc_solvers/c99/d_aux_ip_hard_libstr.c
src/interfaces/c/fortran_order_interface_libstr.c)
SET_TARGET_PROPERTIES(casadi_hpmpc_X64_SSE3_X64_INTEL_CORE PROPERTIES COMPILE_FLAGS "-DLA_BLASFEO -DTARGET_X64_INTEL_CORE -DTARGET_X64_SSE3 -O2 -fPIC -DBLASFEO -m64 -msse3 -DTARGET_X64_SSE3 ${HPMPC_FLAGS}")

target_link_libraries(casadi_hpmpc_X64_SSE3_X64_INTEL_CORE casadi_blasfeo_X64_INTEL_CORE m)
endif()
if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
casadi_external_library(casadi_hpmpc_X64_SSE3_X64_AMD_BULLDOZER src/auxiliary/d_aux_lib4.c
src/auxiliary/d_aux_extern_depend_lib4.c
src/auxiliary/i_aux.c
src/lqcp_solvers/d_back_ric_rec_libstr.c
src/lqcp_solvers/d_part_cond_libstr.c
src/lqcp_solvers/d_tree_back_ric_rec_libstr.c
src/mpc_solvers/d_ip2_res_hard_libstr.c
src/mpc_solvers/d_tree_ip2_res_hard_libstr.c
src/mpc_solvers/d_res_ip_res_hard_libstr.c
src/mpc_solvers/d_tree_res_ip_res_hard_libstr.c
src/mpc_solvers/c99/d_aux_ip_hard_libstr.c
src/interfaces/c/fortran_order_interface_libstr.c)
SET_TARGET_PROPERTIES(casadi_hpmpc_X64_SSE3_X64_AMD_BULLDOZER PROPERTIES COMPILE_FLAGS "-DLA_BLASFEO -DTARGET_X64_AMD_BULLDOZER -DTARGET_X64_SSE3 -O2 -fPIC -DBLASFEO -m64 -msse3 -DTARGET_X64_SSE3 ${HPMPC_FLAGS}")

target_link_libraries(casadi_hpmpc_X64_SSE3_X64_AMD_BULLDOZER casadi_blasfeo_X64_AMD_BULLDOZER m)
endif()
casadi_external_library(casadi_hpmpc_C99_4X4_GENERIC src/auxiliary/d_aux_lib4.c
src/auxiliary/d_aux_extern_depend_lib4.c
src/auxiliary/i_aux.c
src/lqcp_solvers/d_back_ric_rec_libstr.c
src/lqcp_solvers/d_part_cond_libstr.c
src/lqcp_solvers/d_tree_back_ric_rec_libstr.c
src/mpc_solvers/d_ip2_res_hard_libstr.c
src/mpc_solvers/d_tree_ip2_res_hard_libstr.c
src/mpc_solvers/d_res_ip_res_hard_libstr.c
src/mpc_solvers/d_tree_res_ip_res_hard_libstr.c
src/mpc_solvers/c99/d_aux_ip_hard_libstr.c
src/interfaces/c/fortran_order_interface_libstr.c)
SET_TARGET_PROPERTIES(casadi_hpmpc_C99_4X4_GENERIC PROPERTIES COMPILE_FLAGS "-DLA_BLASFEO -DTARGET_GENERIC -DTARGET_C99_4X4 -O2 -fPIC -DBLASFEO -march=native -DTARGET_C99_4X4 ${HPMPC_FLAGS}")

target_link_libraries(casadi_hpmpc_C99_4X4_GENERIC casadi_blasfeo_GENERIC m)
