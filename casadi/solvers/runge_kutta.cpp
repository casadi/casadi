/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


#include "runge_kutta.hpp"

namespace casadi {

  extern "C"
  int CASADI_INTEGRATOR_RK_EXPORT
      casadi_register_integrator_rk(Integrator::Plugin* plugin) {
    plugin->creator = RungeKutta::creator;
    plugin->name = "rk";
    plugin->doc = RungeKutta::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &RungeKutta::options_;
    plugin->deserialize = &RungeKutta::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_INTEGRATOR_RK_EXPORT casadi_load_integrator_rk() {
    Integrator::registerPlugin(casadi_register_integrator_rk);
  }

  RungeKutta::RungeKutta(const std::string& name, const Function& dae, double t0,
      const std::vector<double>& tout)
      : FixedStepIntegrator(name, dae, t0, tout) {
  }

  RungeKutta::~RungeKutta() {
  }

  void RungeKutta::init(const Dict& opts) {
    // Call the base class init
    FixedStepIntegrator::init(opts);

    // Algebraic variables not supported
    casadi_assert(nz_==0 && nrz_==0,
      "Explicit Runge-Kutta integrators do not support algebraic variables");
  }

  void RungeKutta::setup_step() {
    // Continuous-time dynamics, forward problem
    Function f = get_function("dae");

    // Symbolic inputs
    MX t0 = MX::sym("t0", f.sparsity_in(DYN_T));
    MX h = MX::sym("h");
    MX x0 = MX::sym("x0", f.sparsity_in(DYN_X));
    MX p = MX::sym("p", f.sparsity_in(DYN_P));
    MX u = MX::sym("u", f.sparsity_in(DYN_U));

    // Intermediate variables (does not enter in F_, only in G_)
    MX v = MX::sym("v", x0.size1(), x0.size2() * 3);
    std::vector<MX> x = horzsplit(v, x0.size2());
    casadi_assert_dev(x.size() == 3);

    // Definitions of x
    std::vector<MX> x_def(3);

    // Time points
    std::vector<MX> tt(3);

    // Half a step, 6-th of a step
    MX h_half = h / 2, h_sixth = h / 6;

    // Forward integration
    {
      // Arguments when calling f
      std::vector<MX> f_arg(DYN_NUM_IN);
      std::vector<MX> f_res;
      f_arg[DYN_P] = p;
      f_arg[DYN_U] = u;

      // k1
      f_arg[DYN_T] = t0;
      f_arg[DYN_X] = x0;
      f_res = f(f_arg);
      MX k1 = f_res[DYN_ODE];
      MX k1q = f_res[DYN_QUAD];

      // k2
      tt[0] = f_arg[DYN_T] = t0 + h_half;
      x_def[0] = f_arg[DYN_X] = x0 + h_half * k1;
      f_res = f(f_arg);
      MX k2 = f_res[DYN_ODE];
      MX k2q = f_res[DYN_QUAD];

      // k3
      tt[1] = tt[0];
      x_def[1] = f_arg[DYN_X] = x0 + h_half * k2;
      f_res = f(f_arg);
      MX k3 = f_res[DYN_ODE];
      MX k3q = f_res[DYN_QUAD];

      // k4
      tt[2] = f_arg[DYN_T] = t0 + h;
      x_def[2] = f_arg[DYN_X] = x0 + h * k3;
      f_res = f(f_arg);
      MX k4 = f_res[DYN_ODE];
      MX k4q = f_res[DYN_QUAD];

      // Take step
      MX xf = x0 + h_sixth * (k1 + 2*k2 + 2*k3 + k4);
      MX qf = h_sixth * (k1q + 2*k2q + 2*k3q + k4q);

      // Define discrete time dynamics
      f_arg.resize(STEP_NUM_IN);
      f_arg[STEP_T] = t0;
      f_arg[STEP_H] = h;
      f_arg[STEP_X0] = x0;
      f_arg[STEP_V0] = v;
      f_arg[STEP_P] = p;
      f_arg[STEP_U] = u;
      f_res.resize(STEP_NUM_OUT);
      f_res[STEP_XF] = xf;
      f_res[STEP_QF] = qf;
      f_res[STEP_VF] = horzcat(x_def);
      Function F("stepF", f_arg, f_res,
        {"t", "h", "x0", "v0", "p", "u"}, {"xf", "vf", "qf"});
      set_function(F, F.name(), true);
      if (nfwd_ > 0) create_forward("stepF", nfwd_);
    }

    // Backward integration
    if (nadj_ > 0) {
      // Continuous-time dynamics, backward problem
      Function g = get_function("rdae");

      // Symbolic inputs
      MX rx0 = MX::sym("rx0", g.sparsity_in(BDYN_ADJ_ODE));
      MX rp = MX::sym("rp", g.sparsity_in(BDYN_ADJ_QUAD));

      // Intermediate variables (do not enter in G_)
      MX rv = MX::sym("rv", rx0.size1(), 3 * rx0.size2());
      std::vector<MX> rx_def(3);

      // Arguments when calling g
      std::vector<MX> g_arg(BDYN_NUM_IN);
      std::vector<MX> g_res;
      g_arg[BDYN_P] = p;
      g_arg[BDYN_U] = u;
      g_arg[BDYN_ADJ_QUAD] = rp;

      // k1
      g_arg[BDYN_T] = tt[2];
      g_arg[BDYN_X] = x[2];
      g_arg[BDYN_ADJ_ODE] = rx0;
      g_res = g(g_arg);
      MX k1 = g_res[BDYN_ADJ_X];
      MX k1rq = g_res[BDYN_ADJ_P];
      MX k1uq = g_res[BDYN_ADJ_U];

      // k2
      g_arg[BDYN_T] = tt[1];
      g_arg[BDYN_X] = x[1];
      g_arg[BDYN_ADJ_ODE] = rx_def[2] = rx0 + h_half * k1;
      g_res = g(g_arg);
      MX k2 = g_res[BDYN_ADJ_X];
      MX k2rq = g_res[BDYN_ADJ_P];
      MX k2uq = g_res[BDYN_ADJ_U];

      // k3
      g_arg[BDYN_T] = tt[0];
      g_arg[BDYN_X] = x[0];
      g_arg[BDYN_ADJ_ODE] = rx_def[1] = rx0 + h_half * k2;
      g_res = g(g_arg);
      MX k3 = g_res[BDYN_ADJ_X];
      MX k3rq = g_res[BDYN_ADJ_P];
      MX k3uq = g_res[BDYN_ADJ_U];

      // k4
      g_arg[BDYN_T] = t0;
      g_arg[BDYN_X] = x0;
      g_arg[BDYN_ADJ_ODE] = rx_def[0] = rx0 + h * k3;
      g_res = g(g_arg);
      MX k4 = g_res[BDYN_ADJ_X];
      MX k4rq = g_res[BDYN_ADJ_P];
      MX k4uq = g_res[BDYN_ADJ_U];

      // Take step
      MX rxf = rx0 + h_sixth * (k1 + 2*k2 + 2*k3 + k4);
      MX rqf = h_sixth * (k1rq + 2*k2rq + 2*k3rq + k4rq);
      MX uqf = h_sixth * (k1uq + 2*k2uq + 2*k3uq + k4uq);

      // Define discrete time dynamics
      g_arg.resize(BSTEP_NUM_IN);
      g_arg[BSTEP_T] = t0;
      g_arg[BSTEP_H] = h;
      g_arg[BSTEP_X] = x0;
      g_arg[BSTEP_P] = p;
      g_arg[BSTEP_U] = u;
      g_arg[BSTEP_V] = v;
      g_arg[BSTEP_RX0] = rx0;
      g_arg[BSTEP_RP] = rp;
      g_arg[BSTEP_RV0] = rv;
      g_res.resize(BSTEP_NUM_OUT);
      g_res[BSTEP_RXF] = rxf;
      g_res[BSTEP_RVF] = horzcat(rx_def);
      g_res[BSTEP_RQF] = rqf;
      g_res[BSTEP_UQF] = uqf;
      Function G("stepB", g_arg, g_res,
        {"t", "h", "rx0", "rv0", "rp", "x", "v", "p", "u"},
        {"rxf", "rvf", "rqf", "uqf"});
      set_function(G, G.name(), true);
      if (nfwd_ > 0) create_forward("stepB", nfwd_);
    }
  }

  RungeKutta::RungeKutta(DeserializingStream& s) : FixedStepIntegrator(s) {
    s.version("RungeKutta", 2);
  }

  void RungeKutta::serialize_body(SerializingStream &s) const {
    FixedStepIntegrator::serialize_body(s);
    s.version("RungeKutta", 2);
  }

} // namespace casadi
