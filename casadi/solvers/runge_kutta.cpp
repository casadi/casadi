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

    // Half a step, 6-th of a step
    MX h_half = h / 2, h_sixth = h / 6;

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
    f_arg[DYN_T] = t0 + h_half;
    f_arg[DYN_X] = x0 + h_half * k1;
    f_res = f(f_arg);
    MX k2 = f_res[DYN_ODE];
    MX k2q = f_res[DYN_QUAD];

    // k3
    f_arg[DYN_X] = x0 + h_half * k2;
    f_res = f(f_arg);
    MX k3 = f_res[DYN_ODE];
    MX k3q = f_res[DYN_QUAD];

    // k4
    f_arg[DYN_T] = t0 + h;
    f_arg[DYN_X] = x0 + h * k3;
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
    f_arg[STEP_V0] = MX(0, 1);
    f_arg[STEP_P] = p;
    f_arg[STEP_U] = u;
    f_res.resize(STEP_NUM_OUT);
    f_res[STEP_XF] = xf;
    f_res[STEP_QF] = qf;
    f_res[STEP_VF] = MX(0, 1);
    Function F("step", f_arg, f_res,
      {"t", "h", "x0", "v0", "p", "u"}, {"xf", "vf", "qf"});
    set_function(F, F.name(), true);
    if (nfwd_ > 0) create_forward("step", nfwd_);

    // Backward integration
    if (nadj_ > 0) {
      Function adj_F = F.reverse(nadj_);
      set_function(adj_F, adj_F.name(), true);
      if (nfwd_ > 0) {
        create_forward(adj_F.name(), nfwd_);
      }
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
