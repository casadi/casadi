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
    Function f = get_function("dae");

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
    MX t0 = MX::sym("t0", f.sparsity_in(DYN_T));
    MX h = MX::sym("h");
    MX x0 = MX::sym("x0", f.sparsity_in(DYN_X));
    MX p = MX::sym("p", f.sparsity_in(DYN_P));
    MX u = MX::sym("u", f.sparsity_in(DYN_U));

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
    MX qf = MX::zeros(nq1_);

    // End state
    MX xf = D[0] * x0;

    // For all collocation points
    for (casadi_int j = 1; j < deg_ + 1; ++j) {

      // Evaluate the DAE
      std::vector<MX> f_arg(DYN_NUM_IN);
      f_arg[DYN_T] = tt[j];
      f_arg[DYN_P] = p;
      f_arg[DYN_U] = u;
      f_arg[DYN_X] = x[j];
      f_arg[DYN_Z] = z[j];
      std::vector<MX> f_res = f(f_arg);

      // Get an expression for the state derivative at the collocation point
      MX xp_j = C[0][j] * x0;
      for (casadi_int r = 1; r < deg_ + 1; ++r) {
        xp_j += C[r][j] * x[r];
      }

      // Add collocation equation
      eq.push_back(vec(h * f_res[DYN_ODE] - xp_j));

      // Add the algebraic conditions
      eq.push_back(vec(f_res[DYN_ALG]));

      // Add contribution to the final state
      xf += D[j] * x[j];

      // Add contribution to quadratures
      qf += (B[j] * h) * f_res[DYN_QUAD];
    }

    // Form forward discrete time dynamics
    std::vector<MX> F_in(STEP_NUM_IN);
    F_in[STEP_T] = t0;
    F_in[STEP_H] = h;
    F_in[STEP_X0] = x0;
    F_in[STEP_P] = p;
    F_in[STEP_U] = u;
    F_in[STEP_V0] = v;
    std::vector<MX> F_out(STEP_NUM_OUT);
    F_out[STEP_XF] = xf;
    F_out[STEP_VF] = vertcat(eq);
    F_out[STEP_QF] = qf;
    Function F("implicit_step", F_in, F_out,
      {"t", "h", "x0", "v0", "p", "u"}, {"xf", "vf", "qf"});
    set_function(F, F.name(), true);
  }

  void Collocation::reset(IntegratorMemory* mem, bool first_call) const {
    auto m = static_cast<FixedStepMemory*>(mem);

    // Reset the base classes
    ImplicitFixedStepIntegrator::reset(mem, first_call);

    if (first_call) {
      // Initial guess for v (only non-augmented system part)
      double* v = m->v;
      for (casadi_int d = 0; d < deg_; ++d) {
        casadi_copy(m->x, nx1_, v);
        v += nx1_;
        casadi_copy(m->z, nz1_, v);
        v += nz1_;
      }
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
