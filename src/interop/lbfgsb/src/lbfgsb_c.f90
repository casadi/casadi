module lbfgsb_c
   use iso_c_binding
   implicit none (type,external)

   interface
      !> Main function from the original L-BFGS-B Fortran API.
      subroutine setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, &
         wa, iwa, task, iprint, csave, lsave, isave, dsave)
         character(60)    task, csave
         logical          lsave(4)
         integer          n, m, iprint, nbd(n), iwa(3*n), isave(44)
         double precision f, factr, pgtol, x(n), l(n), u(n), g(n), &
            wa(2*m*n + 5*n + 11*m*m + 8*m), dsave(29)
      end subroutine
   end interface

contains

   !> C interface to @ref setulb.
   subroutine alpaqa_setulb_c(n, m, x, l, u, nbd, f, g, factr, pgtol, &
      wa, iwa, task, iprint, csave, lsave, isave, dsave) bind(C, name="alpaqa_setulb_c")
      integer(c_int), intent(in), value :: n, m
      real(c_double), intent(inout) :: x(n)
      real(c_double), intent(in) :: l(n), u(n)
      integer(c_int), intent(in) :: nbd(n)
      real(c_double), intent(inout) :: f
      real(c_double), intent(inout) :: g(n)
      real(c_double), intent(in), value :: factr
      real(c_double), intent(in), value :: pgtol
      real(c_double), intent(out) :: wa(2*m*n + 5*n + 11*m*m + 8*m)
      integer(c_int), intent(out) :: iwa(3*n)
      character(c_char), intent(inout) :: task(60)
      integer(c_int), intent(in), value :: iprint
      character(c_char), intent(inout) :: csave(60)
      logical(c_bool), intent(inout) :: lsave(4)
      integer(c_int), intent(inout) :: isave(44)
      real(c_double), intent(inout) :: dsave(29)

      character(60) :: task_f, csave_f
      logical :: lsave_f(4)
      integer :: i
      do i=1, 60
         task_f(i:i) = task(i)
         csave_f(i:i) = csave(i)
      end do
      do i=1, 4
         lsave_f(i) = lsave(i)
      end do
      call setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, &
         wa, iwa, task_f, iprint, csave_f, lsave_f, isave, dsave)
      do i=1, 4
         lsave(i) = lsave_f(i)
      end do
      do i = 1, 60
         csave(i) = csave_f(i:i)
         task(i) = task_f(i:i)
      end do
   end

end module lbfgsb_c
