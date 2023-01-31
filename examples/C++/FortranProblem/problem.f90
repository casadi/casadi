module problem

   use, intrinsic::iso_c_binding, only: c_double, c_ptrdiff_t
   implicit none (type,external)

   integer(c_ptrdiff_t), parameter :: n = 2
   integer(c_ptrdiff_t), parameter :: m = 1
   real(c_double) :: Q(n, n)
   data Q /3, -1, -1, 3/
   real(c_double) :: A(m, n)
   data A /2, 1/

contains

   ! Evaluate the cost, f(x)
   pure real(c_double) function problem_eval_f(x) bind(C)
      real(c_double), intent(in) :: x(n)

      problem_eval_f = 0.5_c_double * dot_product(x, matmul(Q, x))
   end function

   ! Evaluate the gradient of the cost, ∇f(x)
   pure subroutine problem_eval_grad_f(x, grad_fx) bind(C)
      real(c_double), intent(in) :: x(n)
      real(c_double), intent(out) :: grad_fx(n)

      grad_fx = matmul(Q, x)
   end subroutine

   ! Evaluate the constraints, g(x)
   pure subroutine problem_eval_g(x, gx) bind(C)
      real(c_double), intent(in) :: x(n)
      real(c_double), intent(out) :: gx(m)

      gx = matmul(A, x)
   end subroutine

   ! Evaluate the matrix-vector product of the gradient of the constraints and
   ! the vector y, ∇g(x) y
   pure subroutine problem_eval_grad_g_prod(x, y, grad_gxy) bind(C)
      real(c_double), intent(in) :: x(n)
      real(c_double), intent(in) :: y(m)
      real(c_double), intent(out) :: grad_gxy(n)

      grad_gxy = matmul(transpose(A), y)
   end subroutine

   ! Get the number of variables of the problem
   pure integer(c_ptrdiff_t) function problem_get_num_vars() bind(C)
      problem_get_num_vars = n
   end function

   ! Get the number of constraints of the problem
   pure integer(c_ptrdiff_t) function problem_get_num_constr() bind(C)
      problem_get_num_constr = m
   end function

end module
