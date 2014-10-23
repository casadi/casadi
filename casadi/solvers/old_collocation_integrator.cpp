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


#include "old_collocation_integrator.hpp"
#include "casadi/core/std_vector_tools.hpp"
#include "casadi/core/polynomial.hpp"
#include "casadi/core/matrix/sparsity_tools.hpp"
#include "casadi/core/matrix/matrix_tools.hpp"
#include "casadi/core/sx/sx_tools.hpp"
#include "casadi/core/function/sx_function.hpp"
#include "casadi/core/mx/mx_tools.hpp"

#include "casadi/core/profiling.hpp"
#include "casadi/core/casadi_options.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_INTEGRATOR_OLDCOLLOCATION_EXPORT
  casadi_register_integrator_oldcollocation(IntegratorInternal::Plugin* plugin) {
    plugin->creator = OldCollocationIntegrator::creator;
    plugin->name = "oldcollocation";
    plugin->doc = OldCollocationIntegrator::meta_doc.c_str();
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_INTEGRATOR_OLDCOLLOCATION_EXPORT casadi_load_integrator_oldcollocation() {
    IntegratorInternal::registerPlugin(casadi_register_integrator_oldcollocation);
  }

  OldCollocationIntegrator::OldCollocationIntegrator(const Function& f,
                                                                     const Function& g)
      : IntegratorInternal(f, g) {
    addOption("number_of_finite_elements",     OT_INTEGER,  20,
              "Number of finite elements");
    addOption("interpolation_order",           OT_INTEGER,  3,
              "Order of the interpolating polynomials");
    addOption("collocation_scheme",            OT_STRING,  "radau",
              "Collocation scheme", "radau|legendre");
    addOption("implicit_solver",               OT_STRING,  GenericType(),
              "An implicit function solver");
    addOption("implicit_solver_options",       OT_DICTIONARY, GenericType(),
              "Options to be passed to the implicit solver");
    addOption("expand_f",                      OT_BOOLEAN,  false,
              "Expand the ODE/DAE residual function in an SX graph");
    addOption("expand_q",                      OT_BOOLEAN,  false,
              "Expand the quadrature function in an SX graph");
    addOption("hotstart",                      OT_BOOLEAN,  true,
              "Initialize the trajectory at the previous solution");
    addOption("startup_integrator",            OT_STRING,  GenericType(),
              "An ODE/DAE integrator that can be used to generate a startup trajectory");
    addOption("startup_integrator_options",    OT_DICTIONARY, GenericType(),
              "Options to be passed to the startup integrator");
    setOption("name", "unnamed_old_collocation_integrator");
  }

  void OldCollocationIntegrator::deepCopyMembers(
      std::map<SharedObjectNode*, SharedObject>& already_copied) {
    IntegratorInternal::deepCopyMembers(already_copied);
    startup_integrator_ = deepcopy(startup_integrator_, already_copied);
    implicit_solver_ = deepcopy(implicit_solver_, already_copied);
  }

  OldCollocationIntegrator::~OldCollocationIntegrator() {
  }

  void OldCollocationIntegrator::init() {

    // Call the base class init
    IntegratorInternal::init();

    // Read options
    bool expand_f = getOption("expand_f");

    // Hotstart?
    hotstart_ = getOption("hotstart");

    // Number of finite elements
    int nk = getOption("number_of_finite_elements");

    // Interpolation order
    int deg = getOption("interpolation_order");

    // All collocation time points
    std::vector<double> tau_root = collocationPoints(deg, getOption("collocation_scheme"));

    // Size of the finite elements
    double h = (tf_-t0_)/nk;

    // MX version of the same
    MX h_mx = h;

    // Coefficients of the collocation equation
    vector<vector<MX> > C(deg+1, vector<MX>(deg+1));

    // Coefficients of the collocation equation as DMatrix
    DMatrix C_num = DMatrix::zeros(deg+1, deg+1);

    // Coefficients of the continuity equation
    vector<MX> D(deg+1);

    // Coefficients of the collocation equation as DMatrix
    DMatrix D_num = DMatrix::zeros(deg+1);

    // Coefficients of the quadratures
    DMatrix Q = DMatrix::zeros(deg+1);

    // For all collocation points
    for (int j=0; j<deg+1; ++j) {

      // Construct Lagrange polynomials to get the polynomial basis at the collocation point
      Polynomial p = 1;
      for (int r=0; r<deg+1; ++r) {
        if (r!=j) {
          p *= Polynomial(-tau_root[r], 1)/(tau_root[j]-tau_root[r]);
        }
      }

      // Evaluate the polynomial at the final time to get the coefficients
      // of the continuity equation
      D_num(j) = p(1.0);
      D[j] = D_num(j);

      // Evaluate the time derivative of the polynomial at all collocation points to get
      // the coefficients of the continuity equation
      Polynomial dp = p.derivative();
      for (int r=0; r<deg+1; ++r) {
        C_num(j, r) = dp(tau_root[r]);
        C[j][r] = C_num(j, r);
      }

      // Integrate polynomial to get the coefficients of the quadratures
      Polynomial ip = p.anti_derivative();
      Q(j) = ip(1.0);
    }

    C_num(std::vector<int>(1, 0), ALL) = 0;
    C_num(0, 0)   = 1;

    casadi_assert_message(fabs(sumAll(Q)-1).at(0)<1e-9, "Check on quadrature coefficients");
    casadi_assert_message(fabs(sumAll(D_num)-1).at(0)<1e-9, "Check on collocation coefficients");

    // Initial state
    MX X0 = MX::sym("X0", x0().sparsity());

    // Parameters
    MX P = MX::sym("P", p().sparsity());

    // Backward state
    MX RX0 = MX::sym("RX0", rx0().sparsity());

    // Backward parameters
    MX RP = MX::sym("RP", rp().sparsity());

    // Collocated differential states and algebraic variables
    int nX = (nk*(deg+1)+1)*(nx_+nrx_);
    int nZ = nk*deg*(nz_+nrz_);

    // Unknowns
    MX V = MX::sym("V", nX+nZ);
    int offset = 0;

    // Get collocated states, algebraic variables and times
    vector<vector<MX> > X(nk+1);
    vector<vector<MX> > RX(nk+1);
    vector<vector<MX> > Z(nk);
    vector<vector<MX> > RZ(nk);
    coll_time_.resize(nk+1);
    for (int k=0; k<nk+1; ++k) {
      // Number of time points
      int nj = k==nk ? 1 : deg+1;

      // Allocate differential states expressions at the time points
      X[k].resize(nj);
      RX[k].resize(nj);
      coll_time_[k].resize(nj);

      // Allocate algebraic variable expressions at the collocation points
      if (k!=nk) {
        Z[k].resize(nj-1);
        RZ[k].resize(nj-1);
      }

      // For all time points
      for (int j=0; j<nj; ++j) {
        // Get expressions for the differential state
        X[k][j] = reshape(V(range(offset, offset+nx_)), x0().shape());
        offset += nx_;
        RX[k][j] = reshape(V(range(offset, offset+nrx_)), rx0().shape());
        offset += nrx_;

        // Get the local time
        coll_time_[k][j] = t0_ + h*(k + tau_root[j]);

        // Get expressions for the algebraic variables
        if (j>0) {
          Z[k][j-1] = reshape(V(range(offset, offset+nz_)), z0().shape());
          offset += nz_;
          RZ[k][j-1] = reshape(V(range(offset, offset+nrz_)), rz0().shape());
          offset += nrz_;
        }
      }
    }

    // Check offset for consistency
    casadi_assert(offset==V.size());

    // Constraints
    vector<MX> g;
    g.reserve(2*(nk+1));

    // Quadrature expressions
    MX QF = MX::zeros(qf().sparsity());
    MX RQF = MX::zeros(rqf().sparsity());

    // Counter
    int jk = 0;

    // Add initial condition
    g.push_back(vec(X[0][0]-X0));

    // For all finite elements
    for (int k=0; k<nk; ++k, ++jk) {

      // For all collocation points
      for (int j=1; j<deg+1; ++j, ++jk) {
        // Get the time
        MX tkj = coll_time_[k][j];

        // Get an expression for the state derivative at the collocation point
        MX xp_jk = 0;
        for (int j2=0; j2<deg+1; ++j2) {
          xp_jk += C[j2][j]*X[k][j2];
        }

        // Add collocation equations to the list of equations
        vector<MX> f_in(DAE_NUM_IN);
        f_in[DAE_T] = tkj;
        f_in[DAE_P] = P;
        f_in[DAE_X] = X[k][j];
        f_in[DAE_Z] = Z[k][j-1];

        vector<MX> f_out;
        f_out = f_.call(f_in);
        g.push_back(vec(h_mx*f_out[DAE_ODE] - xp_jk));

        // Add the algebraic conditions
        if (nz_>0) {
          g.push_back(vec(f_out[DAE_ALG]));
        }

        // Add the quadrature
        if (nq_>0) {
          QF += Q(j)*h_mx*f_out[DAE_QUAD];
        }

        // Now for the backward problem
        if (nrx_>0) {

          // Get an expression for the state derivative at the collocation point
          MX rxp_jk = 0;
          for (int j2=0; j2<deg+1; ++j2) {
            rxp_jk += C[j2][j]*RX[k][j2];
          }

          // Add collocation equations to the list of equations
          vector<MX> g_in(RDAE_NUM_IN);
          g_in[RDAE_T] = tkj;
          g_in[RDAE_X] = X[k][j];
          g_in[RDAE_Z] = Z[k][j-1];
          g_in[RDAE_P] = P;
          g_in[RDAE_RP] = RP;
          g_in[RDAE_RX] = RX[k][j];
          g_in[RDAE_RZ] = RZ[k][j-1];

          vector<MX> g_out;
          g_out = g_.call(g_in);
          g.push_back(vec(h_mx*g_out[RDAE_ODE] + rxp_jk));

          // Add the algebraic conditions
          if (nrz_>0) {
            g.push_back(vec(g_out[RDAE_ALG]));
          }

          // Add the backward quadrature
          if (nrq_>0) {
            RQF += Q(j)*h_mx*g_out[RDAE_QUAD];
          }
        }
      }

      // Get an expression for the state at the end of the finite element
      MX xf_k = 0;
      for (int j=0; j<deg+1; ++j) {
        xf_k += D[j]*X[k][j];
      }

      // Add continuity equation to the list of equations
      g.push_back(vec(X[k+1][0] - xf_k));

      if (nrx_>0) {
        // Get an expression for the state at the end of the finite element
        MX rxf_k = 0;
        for (int j=0; j<deg+1; ++j) {
          rxf_k += D[j]*RX[k][j];
        }

        // Add continuity equation to the list of equations
        g.push_back(vec(RX[k+1][0] - rxf_k));
      }
    }

    // Add initial condition for the backward integration
    if (nrx_>0) {
      g.push_back(vec(RX[nk][0]-RX0));
    }

    // Constraint expression
    MX gv = vertcat(g);

    casadi_assert(gv.size2()==1);


    // Make sure that the dimension is consistent with the number of unknowns
    casadi_assert_message(gv.size()==V.size(),
                          "Implicit function unknowns and equations do not match");

    // Implicit function
    vector<MX> ifcn_in(1+INTEGRATOR_NUM_IN);
    ifcn_in[0] = V;
    ifcn_in[1+INTEGRATOR_X0] = X0;
    ifcn_in[1+INTEGRATOR_P] = P;
    ifcn_in[1+INTEGRATOR_RX0] = RX0;
    ifcn_in[1+INTEGRATOR_RP] = RP;
    vector<MX> ifcn_out(1+INTEGRATOR_NUM_OUT);
    ifcn_out[0] = gv;
    ifcn_out[1+INTEGRATOR_XF] = X[nk][0];
    ifcn_out[1+INTEGRATOR_QF] = QF;
    ifcn_out[1+INTEGRATOR_RXF] = RX[0][0];
    ifcn_out[1+INTEGRATOR_RQF] = RQF;
    Function ifcn = MXFunction(ifcn_in, ifcn_out);
    std::stringstream ss_ifcn;
    ss_ifcn << "collocation_implicit_residual_" << getOption("name");
    ifcn.setOption("name", ss_ifcn.str());
    ifcn.init();
    if (expand_f) {
      ifcn = SXFunction(shared_cast<MXFunction>(ifcn));
      ifcn.setOption("name", ss_ifcn.str());
      ifcn.init();
    }

    // Get the root-finding solver creator function
    std::string implicit_function_name = getOption("implicit_solver");

    // Allocate a root-finding solver
    implicit_solver_ = ImplicitFunction(implicit_function_name,
                                        ifcn, Function(), LinearSolver());
    std::stringstream ss_implicit_solver;
    ss_implicit_solver << "collocation_implicitsolver_" << getOption("name");
    implicit_solver_.setOption("name", ss_implicit_solver.str());

    // Pass options
    if (hasSetOption("implicit_solver_options")) {
      const Dictionary& implicit_solver_options = getOption("implicit_solver_options");
      implicit_solver_.setOption(implicit_solver_options);
    }

    // Initialize the solver
    implicit_solver_.init();

    if (hasSetOption("startup_integrator")) {

      // Create the linear solver
      std::string startup_integrator_name = getOption("startup_integrator");

      // Allocate a root-finding solver
      startup_integrator_ = Integrator(startup_integrator_name, f_, g_);

      // Pass options
      startup_integrator_.setOption("t0", coll_time_.front().front());
      startup_integrator_.setOption("tf", coll_time_.back().back());

      std::stringstream ss_startup_integrator;
      ss_startup_integrator << "collocation_startup_" << getOption("name");
      startup_integrator_.setOption("name", ss_startup_integrator.str());

      if (hasSetOption("startup_integrator_options")) {
        const Dictionary& startup_integrator_options = getOption("startup_integrator_options");
        startup_integrator_.setOption(startup_integrator_options);
      }

      // Initialize the startup integrator
      startup_integrator_.init();
    }

    // Mark the system not yet integrated
    integrated_once_ = false;
  }

  void OldCollocationIntegrator::reset() {
    // Set up timers for profiling
    double time_zero=0;
    double time_start=0;
    double time_stop=0;
    if (CasadiOptions::profiling && !CasadiOptions::profilingBinary) {
      time_zero = getRealTime();
      CasadiOptions::profilingLog  << "start " << this << ":" <<getOption("name") << std::endl;
    }

    // Call the base class method
    IntegratorInternal::reset();

    // Pass the inputs
    for (int iind=0; iind<INTEGRATOR_NUM_IN; ++iind) {
      implicit_solver_.input(1+iind).set(input(iind));
    }

    // Pass solution guess (if this is the first integration or if hotstart is disabled)
    if (hotstart_==false || integrated_once_==false) {
      vector<double>& v = implicit_solver_.input().data();

      // Check if an integrator for the startup trajectory has been supplied
      bool has_startup_integrator = !startup_integrator_.isNull();

      // Use supplied integrator, if any
      if (has_startup_integrator) {
        for (int iind=0; iind<INTEGRATOR_NUM_IN; ++iind) {
          startup_integrator_.input(iind).set(input(iind));
        }

        // Reset the integrator
        startup_integrator_.reset();
      }

      // Integrate, stopping at all time points
      int offs=0;
      for (int k=0; k<coll_time_.size(); ++k) {
        for (int j=0; j<coll_time_[k].size(); ++j) {

          if (has_startup_integrator) {
            // Integrate to the time point
            startup_integrator_.integrate(coll_time_[k][j]);
          }

          // Save the differential states
          const DMatrix& x = has_startup_integrator ? startup_integrator_.output(INTEGRATOR_XF) :
              input(INTEGRATOR_X0);
          for (int i=0; i<nx_; ++i) {
            v.at(offs++) = x.at(i);
          }

          // Skip algebraic variables (for now) // FIXME
          if (j>0) {
            if (has_startup_integrator && startup_integrator_.hasOption("init_z")) {
              std::vector<double> init_z = startup_integrator_.getOption("init_z");
              for (int i=0; i<nz_; ++i) {
                v.at(offs++) = init_z.at(i);
              }
            } else {
              offs += nz_;
            }
          }

          // Skip backward states // FIXME
          const DMatrix& rx = input(INTEGRATOR_RX0);
          for (int i=0; i<nrx_; ++i) {
            v.at(offs++) = rx.at(i);
          }

          // Skip backward algebraic variables // FIXME
          if (j>0) {
            offs += nrz_;
          }
        }
      }

      // Print
      if (has_startup_integrator && verbose()) {
        cout << "startup trajectory generated, statistics:" << endl;
        startup_integrator_.printStats();
      }

    }

    if (CasadiOptions::profiling) {
      time_start = getRealTime(); // Start timer
    }

    // Solve the system of equations
    implicit_solver_.evaluate();

    // Save the result
    implicit_solver_.output().set(implicit_solver_.input());

    // Write out profiling information
    if (CasadiOptions::profiling && !CasadiOptions::profilingBinary) {
      time_stop = getRealTime(); // Stop timer
      CasadiOptions::profilingLog
        << (time_stop-time_start)*1e6 << " ns | "
        << (time_stop-time_zero)*1e3 << " ms | "
        << this << ":" << getOption("name")
        << ":0|" << implicit_solver_.get()
        << ":" << implicit_solver_.getOption("name")
        << "|solve system" << std::endl;
    }

    // Mark the system integrated at least once
    integrated_once_ = true;
  }

  void OldCollocationIntegrator::integrate(double t_out) {
    for (int oind=0; oind<INTEGRATOR_NUM_OUT; ++oind) {
      output(oind).set(implicit_solver_.output(1+oind));
    }
  }

  void OldCollocationIntegrator::integrateB(double t_out) {
  }

} // namespace casadi
