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

    // Algebraic variables not supported
    casadi_assert(nz_==0 && nrz_==0,
      "Explicit Runge-Kutta integrators do not support algebraic variables");
  }

  void RungeKutta::setup_step() {
    Function f = create_function(nonaug_oracle_, "odeF",
      {"t", "x", "p", "u"}, {"ode", "quad"});

    // Symbolic inputs
    MX t0 = MX::sym("t0", this->t());
    MX h = MX::sym("h");
    MX x0 = MX::sym("x0", this->x1());
    MX p = MX::sym("p", this->p1());
    MX u = MX::sym("u", this->u1());

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
      std::vector<MX> f_arg(ODE_NUM_IN);
      std::vector<MX> f_res;
      f_arg[ODE_P] = p;
      f_arg[ODE_U] = u;

      // k1
      f_arg[ODE_T] = t0;
      f_arg[ODE_X] = x0;
      f_res = f(f_arg);
      MX k1 = f_res[ODE_ODE];
      MX k1q = f_res[ODE_QUAD];

      // k2
      tt[0] = f_arg[ODE_T] = t0 + h_half;
      x_def[0] = f_arg[ODE_X] = x0 + h_half * k1;
      f_res = f(f_arg);
      MX k2 = f_res[ODE_ODE];
      MX k2q = f_res[ODE_QUAD];

      // k3
      tt[1] = tt[0];
      x_def[1] = f_arg[ODE_X] = x0 + h_half * k2;
      f_res = f(f_arg);
      MX k3 = f_res[ODE_ODE];
      MX k3q = f_res[ODE_QUAD];

      // k4
      tt[2] = f_arg[ODE_T] = t0 + h;
      x_def[2] = f_arg[ODE_X] = x0 + h * k3;
      f_res = f(f_arg);
      MX k4 = f_res[ODE_ODE];
      MX k4q = f_res[ODE_QUAD];

      // Take step
      MX xf = x0 + h_sixth * (k1 + 2*k2 + 2*k3 + k4);
      MX qf = h_sixth * (k1q + 2*k2q + 2*k3q + k4q);

      // Define discrete time dynamics
      f_arg.resize(FSTEP_NUM_IN);
      f_arg[FSTEP_T] = t0;
      f_arg[FSTEP_H] = h;
      f_arg[FSTEP_X0] = x0;
      f_arg[FSTEP_V0] = v;
      f_arg[FSTEP_P] = p;
      f_arg[FSTEP_U] = u;
      f_res.resize(FSTEP_NUM_OUT);
      f_res[FSTEP_XF] = xf;
      f_res[FSTEP_QF] = qf;
      f_res[FSTEP_VF] = horzcat(x_def);
      Function F("stepF", f_arg, f_res,
        {"t", "h", "x0", "v0", "p", "u"}, {"xf", "vf", "qf"});
      set_function(F, F.name(), true);
      if (ns_ > 0) create_forward("stepF", ns_);
    }

    // Backward integration
    if (nrx1_ > 0) {
      // Continuous-time dynamics, backward problem
      Function g = create_function(nonaug_oracle_, "odeB",
        {"t", "x", "p", "u", "rx", "rp"},
        {"rode", "rquad", "uquad"});

      // Symbolic inputs
      MX rx0 = MX::sym("rx0", this->rx1());
      MX rp = MX::sym("rp", this->rp1());

      // Intermediate variables (do not enter in G_)
      MX rv = MX::sym("rv", rx0.size1(), 3 * rx0.size2());
      std::vector<MX> rx_def(3);

      // Arguments when calling g
      std::vector<MX> g_arg(RODE_NUM_IN);
      std::vector<MX> g_res;
      g_arg[RODE_P] = p;
      g_arg[RODE_U] = u;
      g_arg[RODE_RP] = rp;

      // k1
      g_arg[RODE_T] = tt[2];
      g_arg[RODE_X] = x[2];
      g_arg[RODE_RX] = rx0;
      g_res = g(g_arg);
      MX k1 = g_res[RODE_RODE];
      MX k1rq = g_res[RODE_RQUAD];
      MX k1uq = g_res[RODE_UQUAD];

      // k2
      g_arg[RODE_T] = tt[1];
      g_arg[RODE_X] = x[1];
      g_arg[RODE_RX] = rx_def[2] = rx0 + h_half * k1;
      g_res = g(g_arg);
      MX k2 = g_res[RODE_RODE];
      MX k2rq = g_res[RODE_RQUAD];
      MX k2uq = g_res[RODE_UQUAD];

      // k3
      g_arg[RODE_T] = tt[0];
      g_arg[RODE_X] = x[0];
      g_arg[RODE_RX] = rx_def[1] = rx0 + h_half * k2;
      g_res = g(g_arg);
      MX k3 = g_res[RODE_RODE];
      MX k3rq = g_res[RODE_RQUAD];
      MX k3uq = g_res[RODE_UQUAD];

      // k4
      g_arg[RODE_T] = t0;
      g_arg[RODE_X] = x0;
      g_arg[RODE_RX] = rx_def[0] = rx0 + h * k3;
      g_res = g(g_arg);
      MX k4 = g_res[RODE_RODE];
      MX k4rq = g_res[RODE_RQUAD];
      MX k4uq = g_res[RODE_UQUAD];

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
      if (ns_ > 0) create_forward("stepB", ns_);
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
