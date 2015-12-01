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


#include "gurobi_interface.hpp"
#include "casadi/core/function/qpsol.hpp"
#include "casadi/core/std_vector_tools.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_QPSOL_GUROBI_EXPORT
  casadi_register_qpsol_gurobi(Qpsol::Plugin* plugin) {
    plugin->creator = GurobiInterface::creator;
    plugin->name = "gurobi";
    plugin->doc = GurobiInterface::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_QPSOL_GUROBI_EXPORT casadi_load_qpsol_gurobi() {
    Qpsol::registerPlugin(casadi_register_qpsol_gurobi);
  }

  GurobiInterface::GurobiInterface(const std::string& name,
                                   const std::map<std::string, Sparsity>& st)
    : Qpsol(name, st) {

    env_ = 0;
  }

  GurobiInterface::~GurobiInterface() {
    if (env_) {
      GRBfreeenv(env_);
      env_ = 0;
    }
  }

  void GurobiInterface::init() {
    // Initialize the base classes
    Qpsol::init();

    // Load environment
    int flag = GRBloadenv(&env_, 0); // no log file
    casadi_assert_message(!flag && env_, "Failed to create GUROBI environment");

    // Temporary memory
    alloc_w(n_, true); // val
    alloc_iw(n_, true); // ind
    alloc_iw(n_, true); // ind2
  }

  void GurobiInterface::eval(const double** arg, double** res, int* iw, double* w, void* mem) {
    // Inputs
    const double *h=arg[QPSOL_H],
      *g=arg[QPSOL_G],
      *a=arg[QPSOL_A],
      *lba=arg[QPSOL_LBA],
      *uba=arg[QPSOL_UBA],
      *lbx=arg[QPSOL_LBX],
      *ubx=arg[QPSOL_UBX],
      *x0=arg[QPSOL_X0],
      *lam_x0=arg[QPSOL_LAM_X0];

    // Outputs
    double *x=res[QPSOL_X],
      *cost=res[QPSOL_COST],
      *lam_a=res[QPSOL_LAM_A],
      *lam_x=res[QPSOL_LAM_X];

    // Temporary memory
    double *val=w; w+=n_;
    int *ind=iw; iw+=n_;
    int *ind2=iw; iw+=n_;

    // Greate an empty model
    GRBmodel *model = 0;
    try {
      int flag = GRBnewmodel(env_, &model, name_.c_str(), 0, 0, 0, 0, 0, 0);
      casadi_assert_message(!flag, GRBgeterrormsg(env_));

      // Add variables
      for (int i=0; i<n_; ++i) {
        // Get bounds
        double lb = lbx ? lbx[i] : 0., ub = ubx ? ubx[i] : 0.;
        if (isinf(lb)) lb = -GRB_INFINITY;
        if (isinf(ub)) ub =  GRB_INFINITY;

        // Pass to model
        flag = GRBaddvar(model, 0, 0, 0, g ? g[i] : 0., lb, ub, GRB_CONTINUOUS, 0);
        casadi_assert_message(!flag, GRBgeterrormsg(env_));
      }
      flag = GRBupdatemodel(model);
      casadi_assert_message(!flag, GRBgeterrormsg(env_));

      // Add quadratic terms
      const int *H_colind=sparsity_in(QPSOL_H).colind(), *H_row=sparsity_in(QPSOL_H).row();
      for (int i=0; i<n_; ++i) {

        // Quadratic term nonzero indices
        int numqnz = H_colind[1]-H_colind[0];
        casadi_copy(H_row, numqnz, ind);
        H_colind++;
        H_row += numqnz;

        // Corresponding column
        casadi_fill(ind2, numqnz, i);

        // Quadratic term nonzeros
        if (h) {
          casadi_copy(h, numqnz, val);
          casadi_scal(numqnz, 0.5, val);
          h += numqnz;
        } else {
          casadi_fill(val, numqnz, 0.);
        }

        // Pass to model
        flag = GRBaddqpterms(model, numqnz, ind, ind2, val);
        casadi_assert_message(!flag, GRBgeterrormsg(env_));
      }

      // Add constraints
      const int *A_colind=sparsity_in(QPSOL_A).colind(), *A_row=sparsity_in(QPSOL_A).row();
      for (int i=0; i<nc_; ++i) {
        // Get bounds
        double lb = lba ? lba[i] : 0., ub = uba ? uba[i] : 0.;
//        if (isinf(lb)) lb = -GRB_INFINITY;
//        if (isinf(ub)) ub =  GRB_INFINITY;

        // Constraint nonzero indices
        int numnz = A_colind[1]-A_colind[0];
        casadi_copy(A_row, numnz, ind);
        A_colind++;
        A_row += numnz;

        // Constraint nonzeros
        if (a) {
          casadi_copy(a, numnz, val);
          a += numnz;
        } else {
          casadi_fill(val, numnz, 0.);
        }

        // Pass to model
        if (isinf(lb)) {
          if (isinf(ub)) {
            // Neither upper or lower bounds, skip
          } else {
            // Only upper bound
            flag = GRBaddconstr(model, numnz, ind, val, GRB_LESS_EQUAL, ub, 0);
            casadi_assert_message(!flag, GRBgeterrormsg(env_));
          }
        } else {
          if (isinf(ub)) {
            // Only lower bound
            flag = GRBaddconstr(model, numnz, ind, val, GRB_GREATER_EQUAL, lb, 0);
            casadi_assert_message(!flag, GRBgeterrormsg(env_));
          } else if (lb==ub) {
            // Upper and lower bounds equal
            flag = GRBaddconstr(model, numnz, ind, val, GRB_EQUAL, lb, 0);
            casadi_assert_message(!flag, GRBgeterrormsg(env_));
          } else {
            // Both upper and lower bounds
            flag = GRBaddrangeconstr(model, numnz, ind, val, lb, ub, 0);
            casadi_assert_message(!flag, GRBgeterrormsg(env_));
          }
        }
      }

      // Solve the optimization problem
      flag = GRBoptimize(model);
      casadi_assert_message(!flag, GRBgeterrormsg(env_));
      int optimstatus;
      flag = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
      casadi_assert_message(!flag, GRBgeterrormsg(env_));

      // Get the objective value, if requested
      if (cost) {
        flag = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, cost);
        casadi_assert_message(!flag, GRBgeterrormsg(env_));
      }

      // Get the optimal solution, if requested
      if (x) {
        flag = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, n_, x);
        casadi_assert_message(!flag, GRBgeterrormsg(env_));
      }

      // Free memory
      GRBfreemodel(model);

    } catch (...) {
      // Free memory
      if (model) GRBfreemodel(model);
      throw;
    }
  }


} // namespace casadi
