# This is a alternative version of the original CMakeLists.txt
# Joel Andersson
# K.U. Leuven 2012

CMAKE_MINIMUM_REQUIRED(VERSION 2.0)

PROJECT(sundials-2.5mod C)
include_directories(include)

SET(SUNDIALS_INCLUDES 
  include/nvector/nvector_serial.h

  include/sundials/sundials_band.h
  include/sundials/sundials_dense.h
  include/sundials/sundials_direct.h
  include/sundials/sundials_iterative.h
#   include/sundials/sundials_lapack.h
  include/sundials/sundials_math.h
  include/sundials/sundials_nvector.h
  include/sundials/sundials_spbcgs.h
  include/sundials/sundials_spgmr.h
  include/sundials/sundials_sptfqmr.h
  include/sundials/sundials_types.h

  include/cvodes/cvodes_band.h
  include/cvodes/cvodes_bandpre.h
  include/cvodes/cvodes_bbdpre.h
  include/cvodes/cvodes_dense.h
  include/cvodes/cvodes_diag.h
  include/cvodes/cvodes_direct.h
  include/cvodes/cvodes.h
#   include/cvodes/cvodes_lapack.h
  include/cvodes/cvodes_spbcgs.h
  include/cvodes/cvodes_spgmr.h
  include/cvodes/cvodes_spils.h
  include/cvodes/cvodes_sptfqmr.h

  include/idas/idas_band.h
  include/idas/idas_bbdpre.h
  include/idas/idas_dense.h
  include/idas/idas_direct.h
  include/idas/idas.h
#   include/idas/idas_lapack.h
  include/idas/idas_spbcgs.h
  include/idas/idas_spgmr.h
  include/idas/idas_spils.h
  include/idas/idas_sptfqmr.h

  include/kinsol/kinsol_band.h
  include/kinsol/kinsol_bbdpre.h
  include/kinsol/kinsol_dense.h
  include/kinsol/kinsol_direct.h
  include/kinsol/kinsol.h
#   include/kinsol/kinsol_lapack.h
  include/kinsol/kinsol_spbcgs.h
  include/kinsol/kinsol_spgmr.h
  include/kinsol/kinsol_spils.h
  include/kinsol/kinsol_sptfqmr.h
)
  
SET(SUNDIALS_SOURCES
  src/nvec_ser/nvector_serial.c

  src/sundials/sundials_band.c
  src/sundials/sundials_dense.c
  src/sundials/sundials_direct.c
  src/sundials/sundials_iterative.c
  src/sundials/sundials_math.c
  src/sundials/sundials_nvector.c
  src/sundials/sundials_spbcgs.c  
  src/sundials/sundials_spgmr.c
  src/sundials/sundials_sptfqmr.c

  src/cvodes/cvodea.c
  src/cvodes/cvodea_io.c
  src/cvodes/cvodes_band.c
  src/cvodes/cvodes_bandpre.c
  src/cvodes/cvodes_bandpre_impl.h
  src/cvodes/cvodes_bbdpre.c
  src/cvodes/cvodes_bbdpre_impl.h  
  src/cvodes/cvodes.c
  src/cvodes/cvodes_dense.c
  src/cvodes/cvodes_diag.c
  src/cvodes/cvodes_diag_impl.h
  src/cvodes/cvodes_direct.c
  src/cvodes/cvodes_direct_impl.h
  src/cvodes/cvodes_impl.h
  src/cvodes/cvodes_io.c
#   src/cvodes/cvodes_lapack.c
  src/cvodes/cvodes_spbcgs.c
  src/cvodes/cvodes_spgmr.c
  src/cvodes/cvodes_spils.c
  src/cvodes/cvodes_spils_impl.h
  src/cvodes/cvodes_sptfqmr.c

  src/idas/idaa.c
  src/idas/idaa_io.c
  src/idas/idas_band.c
  src/idas/idas_bbdpre.c
  src/idas/idas_bbdpre_impl.h
  src/idas/idas.c
  src/idas/idas_dense.c
  src/idas/idas_direct.c
  src/idas/idas_direct_impl.h
  src/idas/idas_ic.c
  src/idas/idas_impl.h
  src/idas/idas_io.c
#   src/idas/idas_lapack.c
  src/idas/idas_spbcgs.c
  src/idas/idas_spgmr.c
  src/idas/idas_spils.c
  src/idas/idas_spils_impl.h
  src/idas/idas_sptfqmr.c
  
  src/kinsol/kinsol_band.c
  src/kinsol/kinsol_bbdpre.c
  src/kinsol/kinsol_bbdpre_impl.h
  src/kinsol/kinsol.c
  src/kinsol/kinsol_dense.c
  src/kinsol/kinsol_direct.c
  src/kinsol/kinsol_direct_impl.h
  src/kinsol/kinsol_impl.h
  src/kinsol/kinsol_io.c
#   src/kinsol/kinsol_lapack.c
  src/kinsol/kinsol_spbcgs.c
  src/kinsol/kinsol_spgmr.c
  src/kinsol/kinsol_spils.c
  src/kinsol/kinsol_spils_impl.h
  src/kinsol/kinsol_sptfqmr.c
)
  
  
if(ENABLE_STATIC)
add_library(sundials-2.5mod STATIC ${SUNDIALS_INCLUDES} ${SUNDIALS_SOURCES})
endif(ENABLE_STATIC)
if(ENABLE_SHARED)
add_library(sundials-2.5mod SHARED ${SUNDIALS_INCLUDES} ${SUNDIALS_SOURCES})
endif(ENABLE_SHARED)
install(TARGETS sundials-2.5mod
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
