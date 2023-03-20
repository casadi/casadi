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


#include "collocation.hpp"
#include "casadi/core/polynomial.hpp"
#include "casadi/core/casadi_misc.hpp"

namespace casadi {

  extern "C"
  int CASADI_INTEGRATOR_COLLOCATION_EXPORT
      casadi_register_integrator_collocation(Integrator::Plugin* plugin) {
    plugin->creator = Collocation::creator;
    plugin->name = "collocation";
    plugin->doc = Collocation::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &Collocation::options_;
    plugin->deserialize = &Collocation::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_INTEGRATOR_COLLOCATION_EXPORT casadi_load_integrator_collocation() {
    Integrator::registerPlugin(casadi_register_integrator_collocation);
  }

  Collocation::Collocation(const std::string& name, const Function& dae,
      double t0, const std::vector<double>& tout)
      : ImplicitFixedStepIntegrator(name, dae, t0, tout) {
  }

  Collocation::~Collocation() {
  }

  const Options Collocation::options_
  = {{&ImplicitFixedStepIntegrator::options_},
     {{"interpolation_order",
       {OT_INT,
        "Order of the interpolating polynomials"}},
      {"collocation_scheme",
       {OT_STRING,
        "Collocation scheme: radau|legendre"}}
     }
  };

  void Collocation::init(const Dict& opts) {
    // Default options
    deg_ = 3;
    collocation_scheme_ = "radau";

    // Read options
    for (auto&& op : opts) {
      if (op.first=="interpolation_order") {
        deg_ = op.second;
      } else if (op.first=="collocation_scheme") {
        collocation_scheme_ = op.second.to_string();
      }
    }

    // Call the base class init
    ImplicitFixedStepIntegrator::init(opts);
  }

  MX Collocation::algebraic_state_init(const MX& x0, const MX& z0) const {
    MX ret = vertcat(x0, z0);
    return repmat(ret, deg_);
  }
  MX Collocation::algebraic_state_output(const MX& Z) const {
    return Z(Slice(Z.size1()-nz_, Z.size1()));
  }

  void Collocation::setup_step() {
    // Continuous-time dynamics, forward problem
    Function f = get_function("dynF");

    // All collocation time points
    std::vector<double> tau_root = collocation_points(deg_, collocation_scheme_);
    tau_root.insert(tau_root.begin(), 0);

    // Coefficients of the collocation equation
    std::vector<std::vector<double> > C(deg_ + 1, std::vector<double>(deg_ + 1, 0));

    // Coefficients of the continuity equation
    std::vector<double> D(deg_ + 1, 0);

    // Coefficients of the quadratures
    std::vector<double> B(deg_ + 1, 0);

    // For all collocation points
    for (casadi_int j = 0; j < deg_ + 1; ++j) {

      // Construct Lagrange polynomials to get the polynomial basis at the collocation point
      Polynomial p = 1;
      for (casadi_int r = 0; r < deg_+1; ++r) {
        if (r!=j) {
          p *= Polynomial(-tau_root[r], 1) / (tau_root[j] - tau_root[r]);
        }
      }

      // Evaluate the polynomial at the final time to get the
      // coefficients of the continuity equation
      if (collocation_scheme_=="radau") {
        D[j] = j==deg_ ? 1 : 0;
      } else {
        D[j] = p(1.0);
      }

      // Evaluate the time derivative of the polynomial at all collocation points to
      // get the coefficients of the continuity equation
      Polynomial dp = p.derivative();
      for (casadi_int r = 0; r < deg_ + 1; ++r) {
        C[j][r] = dp(tau_root[r]);
      }

      // Integrate polynomial to get the coefficients of the quadratures
      Polynomial ip = p.anti_derivative();
      B[j] = ip(1.0);
    }

    // Symbolic inputs
    MX t0 = MX::sym("t0", this->t());
    MX h = MX::sym("h");
    MX x0 = MX::sym("x0", this->x1());
    MX p = MX::sym("p", this->p1());
    MX u = MX::sym("u", this->u1());

    // Implicitly defined variables (z and x)
    MX v = MX::sym("v", deg_ * (nx1_ + nz1_));
    std::vector<casadi_int> v_offset(1, 0);
    for (casadi_int d = 0; d < deg_; ++d) {
      v_offset.push_back(v_offset.back() + nx1_);
      v_offset.push_back(v_offset.back() + nz1_);
    }
    std::vector<MX> vv = vertsplit(v, v_offset);
    std::vector<MX>::const_iterator vv_it = vv.begin();

    // Collocated states
    std::vector<MX> x(deg_ + 1), z(deg_ + 1);
    for (casadi_int d = 1; d <= deg_; ++d) {
      x[d] = *vv_it++;
      z[d] = *vv_it++;
    }
    casadi_assert_dev(vv_it == vv.end());

    // Collocation time points
    std::vector<MX> tt(deg_ + 1);
    for (casadi_int d = 0; d <= deg_; ++d) {
      tt[d] = t0 + h * tau_root[d];
    }

    // Equations that implicitly define v
    std::vector<MX> eq;

    // Quadratures
    MX qf = MX::zeros(this->q1());

    // End state
    MX xf = D[0] * x0;

    // For all collocation points
    for (casadi_int j = 1; j < deg_ + 1; ++j) {

      // Evaluate the DAE
      std::vector<MX> f_arg(FDYN_NUM_IN);
      f_arg[FDYN_T] = tt[j];
      f_arg[FDYN_P] = p;
      f_arg[FDYN_U] = u;
      f_arg[FDYN_X] = x[j];
      f_arg[FDYN_Z] = z[j];
      std::vector<MX> f_res = f(f_arg);

      // Get an expression for the state derivative at the collocation point
      MX xp_j = C[0][j] * x0;
      for (casadi_int r = 1; r < deg_ + 1; ++r) {
        xp_j += C[r][j] * x[r];
      }

      // Add collocation equation
      eq.push_back(vec(h * f_res[FDYN_ODE] - xp_j));

      // Add the algebraic conditions
      eq.push_back(vec(f_res[FDYN_ALG]));

      // Add contribution to the final state
      xf += D[j] * x[j];

      // Add contribution to quadratures
      qf += (B[j] * h) * f_res[FDYN_QUAD];
    }

    // Form forward discrete time dynamics
    std::vector<MX> F_in(FSTEP_NUM_IN);
    F_in[FSTEP_T] = t0;
    F_in[FSTEP_H] = h;
    F_in[FSTEP_X0] = x0;
    F_in[FSTEP_P] = p;
    F_in[FSTEP_U] = u;
    F_in[FSTEP_V0] = v;
    std::vector<MX> F_out(FSTEP_NUM_OUT);
    F_out[FSTEP_XF] = xf;
    F_out[FSTEP_VF] = vertcat(eq);
    F_out[FSTEP_QF] = qf;
    Function F("implicit_stepF", F_in, F_out,
      {"t", "h", "x0", "v0", "p", "u"}, {"xf", "vf", "qf"});
    set_function(F, F.name(), true);

    // Backwards dynamics
    // NOTE: The following is derived so that it will give the exact adjoint
    // sensitivities whenever g is the reverse mode derivative of f.
    if (nrx1_ > 0) {
      // Continuous-time dynamics, backward problem
      Function g = get_function("dynB");

      // Symbolic inputs
      MX rx0 = MX::sym("rx0", this->rx1());
      MX rp = MX::sym("rp", this->rp1());

      // Implicitly defined variables (rz and rx)
      MX rv = MX::sym("v", deg_ * (nrx1_ + nrz1_));
      std::vector<casadi_int> rv_offset(1, 0);
      for (casadi_int d = 0; d < deg_; ++d) {
        rv_offset.push_back(rv_offset.back() + nrx1_);
        rv_offset.push_back(rv_offset.back() + nrz1_);
      }
      std::vector<MX> rvv = vertsplit(rv, rv_offset);
      std::vector<MX>::const_iterator rvv_it = rvv.begin();

      // Collocated states
      std::vector<MX> rx(deg_ + 1), rz(deg_ + 1);
      for (casadi_int d = 1; d <= deg_; ++d) {
        rx[d] = reshape(*rvv_it++, this->rx1().size());
        rz[d] = reshape(*rvv_it++, this->rz1().size());
      }
      casadi_assert_dev(rvv_it == rvv.end());

      // Equations that implicitly define v
      eq.clear();

      // Quadratures
      MX rqf = MX::zeros(this->rq1());
      MX uqf = MX::zeros(this->uq1());

      // End state
      MX rxf = D[0] * rx0;

      // For all collocation points
      for (casadi_int j = 1; j < deg_ + 1; ++j) {

        // Evaluate the backward DAE
        std::vector<MX> g_arg(BDYN_NUM_IN);
        g_arg[BDYN_T] = tt[j];
        g_arg[BDYN_P] = p;
        g_arg[BDYN_U] = u;
        g_arg[BDYN_X] = x[j];
        g_arg[BDYN_Z] = z[j];
        g_arg[BDYN_RX] = rx[j];
        g_arg[BDYN_RZ] = rz[j];
        g_arg[BDYN_RP] = rp;
        std::vector<MX> g_res = g(g_arg);

        // Get an expression for the state derivative at the collocation point
        MX rxp_j = -D[j] * rx0;
        for (casadi_int r = 1; r < deg_ + 1; ++r) {
          rxp_j += (B[r]*C[j][r]) * rx[r];
        }

        // Add collocation equation
        eq.push_back(vec(h * B[j] * g_res[BDYN_RODE] - rxp_j));

        // Add the algebraic conditions
        eq.push_back(vec(g_res[BDYN_RALG]));

        // Add contribution to the final state
        rxf += -B[j] * C[0][j] * rx[j];

        // Add contribution to quadratures
        rqf += h * B[j] * g_res[BDYN_RQUAD];
        uqf += h * B[j] * g_res[BDYN_UQUAD];
      }

      // Form backward discrete time dynamics
      std::vector<MX> G_in(BSTEP_NUM_IN);
      G_in[BSTEP_T] = t0;
      G_in[BSTEP_H] = h;
      G_in[BSTEP_X] = x0;
      G_in[BSTEP_P] = p;
      G_in[BSTEP_U] = u;
      G_in[BSTEP_V] = v;
      G_in[BSTEP_RX0] = rx0;
      G_in[BSTEP_RP] = rp;
      G_in[BSTEP_RV0] = rv;
      std::vector<MX> G_out(BSTEP_NUM_OUT);
      G_out[BSTEP_RXF] = rxf;
      G_out[BSTEP_RVF] = vertcat(eq);
      G_out[BSTEP_RQF] = rqf;
      G_out[BSTEP_UQF] = uqf;
      Function G("implicit_stepB", G_in, G_out,
        {"t", "h", "rx0", "rv0", "rp", "x", "v", "p", "u"},
        {"rxf", "rvf", "rqf", "uqf"});
      set_function(G, G.name(), true);
    }
  }

  void Collocation::reset(IntegratorMemory* mem,
      const double* x, const double* z, const double* p) const {
    auto m = static_cast<FixedStepMemory*>(mem);

    // Reset the base classes
    ImplicitFixedStepIntegrator::reset(mem, x, z, p);

    // Initial guess for v (only non-augmented system part)
    double* v = m->v;
    for (casadi_int d = 0; d < deg_; ++d) {
      casadi_copy(x, nx1_, v);
      v += nx1_;
      casadi_copy(z, nz1_, v);
      v += nz1_;
    }
  }

  void Collocation::resetB(IntegratorMemory* mem,
      const double* rx, const double* rz, const double* rp) const {
    auto m = static_cast<FixedStepMemory*>(mem);

    // Reset the base classes
    ImplicitFixedStepIntegrator::resetB(mem, rx, rz, rp);

    // Initial guess for rv
    double* rv = m->rv;
    for (casadi_int d = 0; d < deg_; ++d) {
      casadi_copy(rx, nrx_, rv);
      rv += nrx_;
      casadi_copy(rz, nrz_, rv);
      rv += nrz_;
    }
  }

  Collocation::Collocation(DeserializingStream& s) : ImplicitFixedStepIntegrator(s) {
    s.version("Collocation", 2);
    s.unpack("Collocation::deg", deg_);
    s.unpack("Collocation::collocation_scheme", collocation_scheme_);
  }

  void Collocation::serialize_body(SerializingStream &s) const {
    ImplicitFixedStepIntegrator::serialize_body(s);
    s.version("Collocation", 2);
    s.pack("Collocation::deg", deg_);
    s.pack("Collocation::collocation_scheme", collocation_scheme_);
  }

} // namespace casadi
