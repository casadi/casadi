/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
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

    // Algebraic variables not (yet?) supported
    casadi_assert(nz_==0 && nrz_==0,
                          "Explicit Runge-Kutta integrators do not support algebraic variables");
  }

  void RungeKutta::setupFG() {
    f_ = create_function("f", {"x", "z", "p", "t"}, {"ode", "alg", "quad"});
    g_ = create_function("g", {"rx", "rz", "rp", "x", "z", "p", "t"},
                              {"rode", "ralg", "rquad"});

    // Symbolic inputs
    MX x0 = MX::sym("x0", this->x());
    MX p = MX::sym("p", this->p());
    MX t = MX::sym("t", this->t());

    // Intermediate variables (does not enter in F_, only in G_)
    MX v = MX::sym("v", x0.size1(), x0.size2()*3);
    std::vector<MX> x = horzsplit(v, x0.size2());
    casadi_assert_dev(x.size()==3);

    // Definitions of x
    std::vector<MX> x_def(3);

    // Time points
    std::vector<MX> tt(3);

    // Forward integration
    {
      // Arguments when calling f
      std::vector<MX> f_arg(FSTEP_NUM_IN);
      std::vector<MX> f_res;
      f_arg[FSTEP_P] = p;

      // k1
      f_arg[FSTEP_T] = t;
      f_arg[FSTEP_X0] = x0;
      f_res = f_(f_arg);
      MX k1 = f_res[FSTEP_XF];
      MX k1q = f_res[FSTEP_QF];

      // k2
      tt[0] = f_arg[FSTEP_T] = t + h_/2;
      x_def[0] = f_arg[FSTEP_X0] = x0 + (h_/2) * k1;
      f_res = f_(f_arg);
      MX k2 = f_res[FSTEP_XF];
      MX k2q = f_res[FSTEP_QF];

      // k3
      tt[1] = tt[0];
      x_def[1] = f_arg[FSTEP_X0] = x0 + (h_/2) * k2;
      f_res = f_(f_arg);
      MX k3 = f_res[FSTEP_XF];
      MX k3q = f_res[FSTEP_QF];

      // k4
      tt[2] = f_arg[FSTEP_T] = t + h_;
      x_def[2] = f_arg[FSTEP_X0] = x0 + h_ * k3;
      f_res = f_(f_arg);
      MX k4 = f_res[FSTEP_XF];
      MX k4q = f_res[FSTEP_QF];

      // Take step
      MX xf = x0 + (h_/6)*(k1 + 2*k2 + 2*k3 + k4);
      MX qf = (h_/6)*(k1q + 2*k2q + 2*k3q + k4q);

      // Define discrete time dynamics
      // TODO(Joel): Change this and make x_k1, x_k2, x_k3 and x_k4 algebraic outputs
      f_arg[FSTEP_T] = t;
      f_arg[FSTEP_X0] = x0;
      f_arg[FSTEP_P] = p;
      f_arg[FSTEP_Z0] = v;
      f_res[FSTEP_XF] = xf;
      f_res[FSTEP_QF] = qf;
      f_res[FSTEP_RES] = horzcat(x_def);
      F_ = Function("dae", f_arg, f_res);
      alloc(F_);
    }

    // Backward integration
    if (!g_.is_null()) {
      // Symbolic inputs
      MX rx0 = MX::sym("rx0", this->rx());
      MX rp = MX::sym("rp", this->rp());

      // Intermediate variables (do not enter in G_)
      MX rv = MX::sym("rv", rx0.size1(), 3*rx0.size2());
      std::vector<MX> rx_def(3);

      // Arguments when calling g
      std::vector<MX> g_arg(BSTEP_NUM_IN);
      std::vector<MX> g_res;
      g_arg[BSTEP_P] = p;
      g_arg[BSTEP_RP] = rp;

      // k1
      g_arg[BSTEP_T] = tt[2];
      g_arg[BSTEP_X] = x[2];
      g_arg[BSTEP_RX0] = rx0;
      g_res = g_(g_arg);
      MX k1 = g_res[BSTEP_RXF];
      MX k1q = g_res[BSTEP_QF];

      // k2
      g_arg[BSTEP_T] = tt[1];
      g_arg[BSTEP_X] = x[1];
      g_arg[BSTEP_RX0] = rx_def[2] = rx0 + (h_/2) * k1;
      g_res = g_(g_arg);
      MX  k2 = g_res[BSTEP_RXF];
      MX  k2q = g_res[BSTEP_QF];

      // k3
      g_arg[BSTEP_T] = tt[0];
      g_arg[BSTEP_X] = x[0];
      g_arg[BSTEP_RX0] = rx_def[1] = rx0 + (h_/2) * k2;
      g_res = g_(g_arg);
      MX k3 = g_res[BSTEP_RXF];
      MX k3q = g_res[BSTEP_QF];

      // k4
      g_arg[BSTEP_T] = t;
      g_arg[BSTEP_X] = x0;
      g_arg[BSTEP_RX0] = rx_def[0] = rx0 + h_ * k3;
      g_res = g_(g_arg);
      MX k4 = g_res[BSTEP_RXF];
      MX  k4q = g_res[BSTEP_QF];

      // Take step
      MX rxf = rx0 + (h_/6)*(k1 + 2*k2 + 2*k3 + k4);
      MX rqf = (h_/6)*(k1q + 2*k2q + 2*k3q + k4q);

      // Define discrete time dynamics
      g_arg[BSTEP_T] = t;
      g_arg[BSTEP_X] = x0;
      g_arg[BSTEP_P] = p;
      g_arg[BSTEP_Z] = v;
      g_arg[BSTEP_RX0] = rx0;
      g_arg[BSTEP_RP] = rp;
      g_arg[BSTEP_RZ0] = rv;
      g_res[BSTEP_RXF] = rxf;
      g_res[BSTEP_QF] = rqf;
      g_res[BSTEP_RES] = horzcat(rx_def);
      G_ = Function("rdae", g_arg, g_res);
      alloc(G_);
    }
  }

  RungeKutta::RungeKutta(DeserializingStream& s) : FixedStepIntegrator(s) {
    s.version("RungeKutta", 1);
    s.unpack("RungeKutta::f", f_);
    s.unpack("RungeKutta::g", g_);
  }

  void RungeKutta::serialize_body(SerializingStream &s) const {
    FixedStepIntegrator::serialize_body(s);
    s.version("RungeKutta", 1);
    s.pack("RungeKutta::f", f_);
    s.pack("RungeKutta::g", g_);
  }

} // namespace casadi
