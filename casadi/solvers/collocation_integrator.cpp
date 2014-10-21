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


#include "collocation_integrator.hpp"
#include "casadi/core/polynomial.hpp"
#include "casadi/core/std_vector_tools.hpp"
#include "casadi/core/matrix/sparsity_tools.hpp"
#include "casadi/core/matrix/matrix_tools.hpp"
#include "casadi/core/sx/sx_tools.hpp"
#include "casadi/core/function/sx_function.hpp"
#include "casadi/core/mx/mx_tools.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_INTEGRATOR_COLLOCATION_EXPORT
      casadi_register_integrator_collocation(IntegratorInternal::Plugin* plugin) {
    plugin->creator = CollocationIntegrator::creator;
    plugin->name = "collocation";
    plugin->doc = CollocationIntegrator::meta_doc.c_str();
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_INTEGRATOR_COLLOCATION_EXPORT casadi_load_integrator_collocation() {
    IntegratorInternal::registerPlugin(casadi_register_integrator_collocation);
  }

  CollocationIntegrator::CollocationIntegrator(const Function& f,
                                                               const Function& g)
      : ImplicitFixedStepIntegrator(f, g) {
    addOption("interpolation_order",           OT_INTEGER,  3,
              "Order of the interpolating polynomials");
    addOption("collocation_scheme",            OT_STRING,  "radau",
              "Collocation scheme", "radau|legendre");
    setOption("name", "unnamed_collocation_integrator");
  }

  void CollocationIntegrator::deepCopyMembers(
      std::map<SharedObjectNode*, SharedObject>& already_copied) {
    ImplicitFixedStepIntegrator::deepCopyMembers(already_copied);
  }

  CollocationIntegrator::~CollocationIntegrator() {
  }

  void CollocationIntegrator::init() {

    // Call the base class init
    ImplicitFixedStepIntegrator::init();

  }

  void CollocationIntegrator::setupFG() {

    // Interpolation order
    deg_ = getOption("interpolation_order");

    // All collocation time points
    std::vector<long double> tau_root = collocationPointsL(deg_, getOption("collocation_scheme"));

    // Coefficients of the collocation equation
    vector<vector<double> > C(deg_+1, vector<double>(deg_+1, 0));

    // Coefficients of the continuity equation
    vector<double> D(deg_+1, 0);

    // Coefficients of the quadratures
    vector<double> B(deg_+1, 0);

    // For all collocation points
    for (int j=0; j<deg_+1; ++j) {

      // Construct Lagrange polynomials to get the polynomial basis at the collocation point
      Polynomial p = 1;
      for (int r=0; r<deg_+1; ++r) {
        if (r!=j) {
          p *= Polynomial(-tau_root[r], 1)/(tau_root[j]-tau_root[r]);
        }
      }

      // Evaluate the polynomial at the final time to get the
      // coefficients of the continuity equation
      D[j] = zeroIfSmall(p(1.0L));

      // Evaluate the time derivative of the polynomial at all collocation points to
      // get the coefficients of the continuity equation
      Polynomial dp = p.derivative();
      for (int r=0; r<deg_+1; ++r) {
        C[j][r] = zeroIfSmall(dp(tau_root[r]));
      }

      // Integrate polynomial to get the coefficients of the quadratures
      Polynomial ip = p.anti_derivative();
      B[j] = zeroIfSmall(ip(1.0L));
    }

    // Symbolic inputs
    MX x0 = MX::sym("x0", f_.input(DAE_X).sparsity());
    MX p = MX::sym("p", f_.input(DAE_P).sparsity());
    MX t = MX::sym("t", f_.input(DAE_T).sparsity());

    // Implicitly defined variables (z and x)
    MX v = MX::sym("v", deg_*(nx_+nz_));
    vector<int> v_offset(1, 0);
    for (int d=0; d<deg_; ++d) {
      v_offset.push_back(v_offset.back()+nx_);
      v_offset.push_back(v_offset.back()+nz_);
    }
    vector<MX> vv = vertsplit(v, v_offset);
    vector<MX>::const_iterator vv_it = vv.begin();

    // Collocated states
    vector<MX> x(deg_+1), z(deg_+1);
    for (int d=1; d<=deg_; ++d) {
      x[d] = reshape(*vv_it++, this->x0().shape());
      z[d] = reshape(*vv_it++, this->z0().shape());
    }
    casadi_assert(vv_it==vv.end());

    // Collocation time points
    vector<MX> tt(deg_+1);
    for (int d=0; d<=deg_; ++d) {
      tt[d] = t + h_*tau_root[d];
    }

    // Equations that implicitly define v
    vector<MX> eq;

    // Quadratures
    MX qf = MX::zeros(f_.output(DAE_QUAD).sparsity());

    // End state
    MX xf = D[0]*x0;

    // For all collocation points
    for (int j=1; j<deg_+1; ++j) {
      //for (int j=deg_; j>=1; --j) {

      // Evaluate the DAE
      vector<MX> f_arg(DAE_NUM_IN);
      f_arg[DAE_T] = tt[j];
      f_arg[DAE_P] = p;
      f_arg[DAE_X] = x[j];
      f_arg[DAE_Z] = z[j];
      vector<MX> f_res = f_.call(f_arg);

      // Get an expression for the state derivative at the collocation point
      MX xp_j = C[0][j] * x0;
      for (int r=1; r<deg_+1; ++r) {
        xp_j += C[r][j] * x[r];
      }

      // Add collocation equation
      eq.push_back(vec(h_*f_res[DAE_ODE] - xp_j));

      // Add the algebraic conditions
      eq.push_back(vec(f_res[DAE_ALG]));

      // Add contribution to the final state
      xf += D[j]*x[j];

      // Add contribution to quadratures
      qf += (B[j]*h_)*f_res[DAE_QUAD];
    }

    // Form forward discrete time dynamics
    vector<MX> F_in(DAE_NUM_IN);
    F_in[DAE_T] = t;
    F_in[DAE_X] = x0;
    F_in[DAE_P] = p;
    F_in[DAE_Z] = v;
    vector<MX> F_out(DAE_NUM_OUT);
    F_out[DAE_ODE] = xf;
    F_out[DAE_ALG] = vertcat(eq);
    F_out[DAE_QUAD] = qf;
    F_ = MXFunction(F_in, F_out);
    F_.init();

    // Backwards dynamics
    // NOTE: The following is derived so that it will give the exact adjoint
    // sensitivities whenever g is the reverse mode derivative of f.
    if (!g_.isNull()) {

      // Symbolic inputs
      MX rx0 = MX::sym("x0", g_.input(RDAE_RX).sparsity());
      MX rp = MX::sym("p", g_.input(RDAE_RP).sparsity());

      // Implicitly defined variables (rz and rx)
      MX rv = MX::sym("v", deg_*(nrx_+nrz_));
      vector<int> rv_offset(1, 0);
      for (int d=0; d<deg_; ++d) {
        rv_offset.push_back(rv_offset.back()+nrx_);
        rv_offset.push_back(rv_offset.back()+nrz_);
      }
      vector<MX> rvv = vertsplit(rv, rv_offset);
      vector<MX>::const_iterator rvv_it = rvv.begin();

      // Collocated states
      vector<MX> rx(deg_+1), rz(deg_+1);
      for (int d=1; d<=deg_; ++d) {
        rx[d] = reshape(*rvv_it++, this->rx0().shape());
        rz[d] = reshape(*rvv_it++, this->rz0().shape());
      }
      casadi_assert(rvv_it==rvv.end());

      // Equations that implicitly define v
      eq.clear();

      // Quadratures
      MX rqf = MX::zeros(g_.output(RDAE_QUAD).sparsity());

      // End state
      MX rxf = D[0]*rx0;

      // For all collocation points
      for (int j=1; j<deg_+1; ++j) {

        // Evaluate the backward DAE
        vector<MX> g_arg(RDAE_NUM_IN);
        g_arg[RDAE_T] = tt[j];
        g_arg[RDAE_P] = p;
        g_arg[RDAE_X] = x[j];
        g_arg[RDAE_Z] = z[j];
        g_arg[RDAE_RX] = rx[j];
        g_arg[RDAE_RZ] = rz[j];
        g_arg[RDAE_RP] = rp;
        vector<MX> g_res = g_.call(g_arg);

        // Get an expression for the state derivative at the collocation point
        MX rxp_j = -D[j]*rx0;
        for (int r=1; r<deg_+1; ++r) {
          rxp_j += (B[r]*C[j][r]) * rx[r];
        }

        // Add collocation equation
        eq.push_back(vec(h_*B[j]*g_res[RDAE_ODE] - rxp_j));

        // Add the algebraic conditions
        eq.push_back(vec(g_res[RDAE_ALG]));

        // Add contribution to the final state
        rxf += -B[j]*C[0][j]*rx[j];

        // Add contribution to quadratures
        rqf += h_*B[j]*g_res[RDAE_QUAD];
      }

      // Form backward discrete time dynamics
      vector<MX> G_in(RDAE_NUM_IN);
      G_in[RDAE_T] = t;
      G_in[RDAE_X] = x0;
      G_in[RDAE_P] = p;
      G_in[RDAE_Z] = v;
      G_in[RDAE_RX] = rx0;
      G_in[RDAE_RP] = rp;
      G_in[RDAE_RZ] = rv;
      vector<MX> G_out(RDAE_NUM_OUT);
      G_out[RDAE_ODE] = rxf;
      G_out[RDAE_ALG] = vertcat(eq);
      G_out[RDAE_QUAD] = rqf;
      G_ = MXFunction(G_in, G_out);
      G_.init();
    }
  }


  double CollocationIntegrator::zeroIfSmall(double x) {
    return fabs(x) < numeric_limits<double>::epsilon() ? 0 : x;
  }

  void CollocationIntegrator::calculateInitialConditions() {
    vector<double>::const_iterator x0_it = input(INTEGRATOR_X0).begin();
    vector<double>::const_iterator z_it = input(INTEGRATOR_Z0).begin();
    vector<double>::iterator Z_it = Z_.begin();
    for (int d=0; d<deg_; ++d) {
      copy(x0_it, x0_it+nx_, Z_it);
      Z_it += nx_;
      copy(z_it, z_it+nz_, Z_it);
      Z_it += nz_;
    }
    casadi_assert(Z_it==Z_.end());
  }

  void CollocationIntegrator::calculateInitialConditionsB() {
    vector<double>::const_iterator rx0_it = input(INTEGRATOR_RX0).begin();
    vector<double>::const_iterator rz_it = input(INTEGRATOR_RZ0).begin();
    vector<double>::iterator RZ_it = RZ_.begin();
    for (int d=0; d<deg_; ++d) {
      copy(rx0_it, rx0_it+nrx_, RZ_it);
      RZ_it += nrx_;
      copy(rz_it, rz_it+nrz_, RZ_it);
      RZ_it += nrz_;
    }
    casadi_assert(RZ_it==RZ_.end());
  }

} // namespace casadi
