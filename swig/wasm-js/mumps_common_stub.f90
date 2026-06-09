! Storage for the MUMPS libseq (sequential fake-MPI) COMMON block.
!
! flang does not emit storage for a Fortran COMMON block in a .o that
! merely references it (unlike gfortran's `.comm` tentative-definition
! convention), so wasm-ld treats `mpif_libseq_` as a hard undefined
! symbol when linking the IPOPT/MUMPS side module.  This BLOCK DATA gives
! the COMMON block exactly one definition.  Compiled in-tree (see
! swig/wasm-js/CMakeLists.txt) so the wasm toolchain image stays generic.
block data mumps_common_storage
  integer :: mpi_in_place
  common /mpif_libseq/ mpi_in_place
  data mpi_in_place /-1/
end block data
