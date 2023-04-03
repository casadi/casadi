//
//    MIT No Attribution
//
//    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
//
//    Permission is hereby granted, free of charge, to any person obtaining a copy of this
//    software and associated documentation files (the "Software"), to deal in the Software
//    without restriction, including without limitation the rights to use, copy, modify,
//    merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
//    permit persons to whom the Software is furnished to do so.
//
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
//    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
//    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
//    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//

// SYMBOL "newton_mem"

template<typename T1>
struct casadi_newton_mem {
  casadi_int n;
  T1 abstol;
  T1 abstol_step;
  T1* x;
  T1* g;
  T1* jac_g_x;

  const casadi_int* sp_a;
  const casadi_int* sp_v;
  const casadi_int* sp_r;
  const casadi_int* prinv;
  const casadi_int* pc;

  T1* lin_w;
  T1* lin_v;
  T1* lin_r;
  T1* lin_beta;
};

// C-REPLACE "casadi_newton_mem<T1>" "struct casadi_newton_mem"
// SYMBOL "newton"
template<typename T1>
int casadi_newton(const casadi_newton_mem<T1>* m) {
    // Check tolerance on residual
    if (m->abstol>0 && casadi_norm_inf(m->n, m->g) <= m->abstol) return 1;

    // Factorize J
    casadi_qr(m->sp_a, m->jac_g_x, m->lin_w,
              m->sp_v,  m->lin_v, m->sp_r, m->lin_r, m->lin_beta,
              m->prinv, m->pc);
    // Solve J^(-1) g
    casadi_qr_solve(m->g, 1, 0, m->sp_v, m->lin_v, m->sp_r, m->lin_r, m->lin_beta,
                    m->prinv, m->pc, m->lin_w);

    // Update Xk+1 = Xk - J^(-1) g
    casadi_axpy(m->n, -1., m->g, m->x);

    // Check tolerance on step
    if (m->abstol_step>0 && casadi_norm_inf(m->n, m->g) <= m->abstol_step) return 2;

    // We will need another newton step
    return 0;
}
