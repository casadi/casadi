! Minimal omp_lib stub for wasm32-emscripten flang build.
! Provides the few intrinsics MUMPS / IPOPT / etc. reference at the
! Fortran level.  Single-threaded semantics: thread_num=0, num=1,
! in_parallel=false, no locking.
module omp_lib
  implicit none
  integer, parameter :: omp_lock_kind = 8
  integer, parameter :: omp_nest_lock_kind = 8
  integer, parameter :: omp_sched_kind = 4
  integer, parameter :: omp_proc_bind_kind = 4
contains
  function omp_get_thread_num() result(n)
    integer :: n
    n = 0
  end function
  function omp_get_num_threads() result(n)
    integer :: n
    n = 1
  end function
  function omp_get_max_threads() result(n)
    integer :: n
    n = 1
  end function
  function omp_in_parallel() result(b)
    logical :: b
    b = .false.
  end function
  function omp_get_wtime() result(t)
    double precision :: t
    t = 0d0
  end function
  subroutine omp_set_num_threads(n)
    integer, intent(in) :: n
  end subroutine
end module
