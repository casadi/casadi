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

// Top level of the acados Integrator plugin, shared by the VM and generated code:
// one evaluation of an AcadosFunction (integrator, jacobian, reverse or jacobian of reverse).
// Inputs and outputs follow the Integrator I/O scheme; requires acados_sim.hpp, acados_chain.hpp

// C-REPLACE "casadi_acados_prob<T1>" "struct casadi_acados_prob"
// C-REPLACE "casadi_acados_data<T1>" "struct casadi_acados_data"
// C-REPLACE "casadi_acados_sim_prob<T1>" "struct casadi_acados_sim_prob"
// C-REPLACE "casadi_acados_sim_data<T1>" "struct casadi_acados_sim_data"
// C-REPLACE "casadi_acados_chain_prob<T1>" "struct casadi_acados_chain_prob"
// C-REPLACE "casadi_acados_chain_data<T1>" "struct casadi_acados_chain_data"
// C-REPLACE "casadi_acados_sim_solve<T1>" "casadi_acados_sim_solve"

template<typename T1>
struct casadi_acados_prob {
  // acados sim and chain over the output intervals (sim.nx, sim.nu set by setup)
  casadi_acados_sim_prob<T1> sim;
  casadi_acados_chain_prob<T1> chain;
  // 0: integrator, 1: jacobian, 2: reverse (nadj directions), 3: jacobian of reverse
  int mode;
  casadi_int nadj;
  // Value of zf: not supported (acados only reports z at the start of an interval)
  T1 nan;
};

template<typename T1>
struct casadi_acados_data {
  const casadi_acados_prob<T1>* prob;
  casadi_acados_sim_data<T1> sim;
  casadi_acados_chain_data<T1> chain;
};

// SYMBOL "acados_setup"
template<typename T1>
void casadi_acados_setup(casadi_acados_prob<T1>* p) {
  // acados controls: [u; p_d]
  p->sim.nx = p->chain.nx;
  p->sim.nu = p->chain.nu + p->chain.npd;
  // Sensitivities needed by the mode
  p->sim.sens_forw = p->mode == 1 || p->mode == 3;
  p->sim.sens_adj = p->mode >= 2;
  p->sim.sens_hess = p->mode == 3;
  p->chain.jac = p->sim.sens_forw;
  p->chain.hess = p->sim.sens_hess;
  casadi_acados_sim_setup(&p->sim);
  casadi_acados_chain_setup(&p->chain);
}

// SYMBOL "acados_work"
template<typename T1>
void casadi_acados_work(const casadi_acados_prob<T1>* p, casadi_int* sz_w) {
  casadi_acados_sim_work(&p->sim, sz_w);
  casadi_acados_chain_work(&p->chain, sz_w);
}

// SYMBOL "acados_set_work"
template<typename T1>
void casadi_acados_set_work(casadi_acados_data<T1>* d, const T1*** arg, T1*** res,
    casadi_int** iw, T1** w) {
  d->sim.prob = &d->prob->sim;
  casadi_acados_sim_set_work(&d->sim, arg, res, iw, w);
  d->chain.prob = &d->prob->chain;
  d->chain.solve = casadi_acados_sim_solve<T1>;
  d->chain.ctx = &d->sim;
  casadi_acados_chain_set_work(&d->chain, arg, res, iw, w);
}

// SYMBOL "acados_eval"
// Returns 0 on success, 1 if acados or a model function failed, 2 if the work is too small
template<typename T1>
int casadi_acados_eval(casadi_acados_data<T1>* d, const T1** arg, T1** res) {
  const casadi_acados_prob<T1>* p = d->prob;
  casadi_acados_chain_data<T1>* c = &d->chain;
  // Integrator I/O: inputs x0, z0, p, u, adj_xf, ...; outputs xf, zf, qf, adj_x0, ..., adj_u
  const casadi_int ni = 7, no = 7, X0 = 0, Z0 = 1, P = 2, U = 3, XF = 0, ZF = 1;
  // Reverse mode: seed on xf is input ni + no + XF; jacobian of reverse has ni + 2*no inputs
  const casadi_int seed_xf = ni + no + XF, nj = ni + 2 * no;
  casadi_int nx = p->chain.nx, nu = p->chain.nu, np = p->chain.np, nt = p->chain.nt;
  casadi_int nr = nx * nt, nW = p->chain.nW, a, k, j;
  // Blocks of w = [x0; u; p_d]: input index, offset, length
  casadi_int ind[3], off[3], len[3];
  const T1 *x0, *pp, *u, *seed;
  T1* r;
  int flag = 0;
  ind[0] = X0; off[0] = 0; len[0] = nx;
  ind[1] = U; off[1] = nx; len[1] = nu * nt;
  ind[2] = P; off[2] = nx + nu * nt; len[2] = p->chain.npd;
  x0 = arg[X0] ? arg[X0] : c->zero;
  pp = arg[P] ? arg[P] : c->zero;
  u = arg[U] ? arg[U] : c->zero;
  d->sim.z0 = arg[Z0];
  if (casadi_acados_sim_init(&d->sim) < 0) return 2;
  if (p->mode == 0) {
    casadi_fill(res[ZF], p->sim.nz * nt, p->nan);
    flag = casadi_acados_chain_nom(c, x0, pp, u, res[XF]);
  } else if (p->mode == 1) {
    flag = casadi_acados_chain_jac(c, x0, pp, u);
    for (k = 0; k < 3; ++k) casadi_copy(c->J + nr * off[k], nr * len[k], res[XF * ni + ind[k]]);
  } else if (p->mode == 2) {
    seed = arg[seed_xf];
    flag = casadi_acados_chain_fwd(c, x0, pp, u, nt - 1, 0);
    for (a = 0; a < p->nadj && !flag; ++a) {
      flag = casadi_acados_chain_adj(c, pp, u, seed ? seed + a * nr : 0,
        res[X0] ? res[X0] + a * nx : 0, res[U] ? res[U] + a * nu * nt : 0,
        res[P] && len[2] ? res[P] + a * np : 0);
    }
  } else {
    flag = casadi_acados_chain_hess(c, x0, pp, u, arg[seed_xf]);
    for (k = 0; k < 3; ++k) {
      for (j = 0; j < 3; ++j) {
        r = res[ind[k] * nj + ind[j]];
        if (r) casadi_acados_block(c->H, nW, off[k], off[j], len[k], len[j], 0, r);
      }
      // d(adj_o)/d(adj_xf) = (d xf / d o)^T
      r = res[ind[k] * nj + seed_xf];
      if (r) casadi_acados_block(c->J, nr, off[k], 0, len[k], nr, 1, r);
    }
  }
  return flag ? 1 : 0;
}
