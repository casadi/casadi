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

#include "cplex_interface.hpp"
#include "casadi/core/std_vector_tools.hpp"
#include "casadi/core/matrix/matrix_tools.hpp"
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <vector>

#include "ilcplex/cplex.h"

namespace casadi {

  using namespace std;

  extern "C"
  int CASADI_QPSOLVER_CPLEX_EXPORT
  casadi_register_qpsolver_cplex(QpSolverInternal::Plugin* plugin) {
    plugin->creator = CplexInterface::creator;
    plugin->name = "cplex";
    plugin->doc = CplexInterface::meta_doc.c_str();
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_QPSOLVER_CPLEX_EXPORT casadi_load_qpsolver_cplex() {
    QpSolverInternal::registerPlugin(casadi_register_qpsolver_cplex);
  }

  CplexInterface::CplexInterface(const std::vector<Sparsity>& st) : QpSolverInternal(st) {
    // Options available
    addOption("qp_method",    OT_STRING, "automatic", "Determines which CPLEX algorithm to use.",
              "automatic|primal_simplex|dual_simplex|network|barrier|sifting|concurrent|crossover");
    addOption("dump_to_file",   OT_BOOLEAN,        false, "Dumps QP to file in CPLEX format.");
    addOption("dump_filename",   OT_STRING,     "qp.dat", "The filename to dump to.");
    addOption("tol",               OT_REAL,         1E-6, "Tolerance of solver");
    addOption("dep_check",      OT_STRING,         "off", "Detect redundant constraints.",
              "automatic:-1|off:0|begin:1|end:2|both:3");
    addOption("simplex_maxiter", OT_INTEGER,   2100000000, "Maximum number of simplex iterations.");
    addOption("barrier_maxiter", OT_INTEGER,   2100000000, "Maximum number of barrier iterations.");
    addOption("warm_start",      OT_BOOLEAN,        false,
              "Use warm start with simplex methods (affects only the simplex methods).");
    addOption("convex",          OT_BOOLEAN,         true,
              "Indicates if the QP is convex or not (affects only the barrier method).");

    // Setting warm-start flag
    is_warm_ = false;

    // Set pointer to zero to avoid deleting a nonexisting instance
    env_ = 0;
    lp_ = 0;
  }
  void CplexInterface::init() {
    // Free any existing Cplex instance
    freeCplex();

    // Call the init method of the base class
    QpSolverInternal::init();

    qp_method_     = getOptionEnumValue("qp_method");
    dump_to_file_  = getOption("dump_to_file");
    tol_           = getOption("tol");
    //  dump_filename_ = getOption("dump_filename");

    int status;
    casadi_assert(env_==0);
    env_ = CPXopenCPLEX(&status);
    casadi_assert_message(env_!=0, "CPLEX: Cannot initialize CPLEX environment. STATUS: "
                          << status);

    // Turn on some debug messages if requested
    if (verbose()) {
      CPXsetintparam(env_, CPX_PARAM_SCRIND, CPX_ON);
    } else {
      CPXsetintparam(env_, CPX_PARAM_SCRIND, CPX_OFF);
    }
    if (status) {
      std::cout << "CPLEX: Problem with setting parameter... ERROR: " << status << std::endl;
    }

    /* SETTING OPTIONS */
    // Optimality tolerance
    status = CPXsetdblparam(env_, CPX_PARAM_EPOPT, tol_);
    // Feasibility tolerance
    status = CPXsetdblparam(env_, CPX_PARAM_EPRHS, tol_);
    // We start with barrier if crossover was chosen.
    if (qp_method_ == 7) {
      status = CPXsetintparam(env_, CPX_PARAM_QPMETHOD, 4);
      // Warm-start is default with this algorithm
      setOption("warm_start", true);
    } else { // Otherwise we just chose the algorithm
      status = CPXsetintparam(env_, CPX_PARAM_QPMETHOD, qp_method_);
    }
    // Setting dependency check option
    status = CPXsetintparam(env_, CPX_PARAM_DEPIND, getOptionEnumValue("dep_check"));
    // Setting barrier iteration limit
    status = CPXsetintparam(env_, CPX_PARAM_BARITLIM, getOption("barrier_maxiter"));
    // Setting simplex iteration limit
    status = CPXsetintparam(env_, CPX_PARAM_ITLIM, getOption("simplex_maxiter"));
    if (qp_method_ == 7) {
      // Setting crossover algorithm
      status = CPXsetintparam(env_, CPX_PARAM_BARCROSSALG, 1);
    }
    if (!static_cast<bool>(getOption("convex"))) {
      // Enabling non-convex QPs
      status = CPXsetintparam(env_, CPX_PARAM_SOLUTIONTARGET, CPX_SOLUTIONTARGET_FIRSTORDER);
    }

    // Exotic parameters, once they might become options...

    // Do careful numerics with numerically unstable problem
    //status = CPXsetintparam(env_, CPX_PARAM_NUMERICALEMPHASIS, 1);
    // Set scaling approach
    //status = CPXsetintparam(env_, CPX_PARAM_SCAIND, 1);
    // Set Markowitz tolerance
    //status = CPXsetdblparam(env_, CPX_PARAM_EPMRK, 0.9);

    // Doing allocation of CPLEX data
    // Objective is to be minimized
    objsen_ = CPX_MIN;

    // Allocation of data
    // Type of constraint
    sense_.resize(nc_);
    // Right-hand side of constraints
    rhs_.resize(nc_);
    // Range value for lower AND  upper bounded constraints
    rngval_.resize(nc_);
    // Basis for primal variables
    cstat_.resize(n_);
    rstat_.resize(nc_);

    // Matrix A, count the number of elements per column
    const Sparsity& A_sp = input(QP_SOLVER_A).sparsity();
    matcnt_.resize(A_sp.size2());
    transform(A_sp.colind().begin()+1, A_sp.colind().end(), A_sp.colind().begin(), matcnt_.begin(),
              minus<int>());

    // Matrix H, count the number of elements per column
    const Sparsity& H_sp = input(QP_SOLVER_H).sparsity();
    qmatcnt_.resize(H_sp.size2());
    transform(H_sp.colind().begin()+1, H_sp.colind().end(), H_sp.colind().begin(), qmatcnt_.begin(),
              minus<int>());

    casadi_assert(lp_==0);
    lp_ = CPXcreateprob(env_, &status, "QP from CasADi");
  }

  void CplexInterface::evaluate() {

    if (inputs_check_) checkInputs();

    int status;

    // We change method in crossover
    if (is_warm_ && qp_method_ == 7) {
      status = CPXsetintparam(env_, CPX_PARAM_QPMETHOD, 1);
    }

    // Looping over constraints
    const vector<double>& lba = input(QP_SOLVER_LBA).data();
    const vector<double>& uba = input(QP_SOLVER_UBA).data();

    for (int i = 0; i < nc_; ++i) {
      // CPX_INFBOUND

      // Equality
      if (uba[i] - lba[i] < 1e-20) {
        sense_[i] = 'E';
        rhs_[i] = lba[i];
        rngval_[i] = 0.;
      } else if (lba[i] < -CPX_INFBOUND) {
        // Ineq - no lower bound
        sense_[i] = 'L';
        rhs_[i] = uba[i];
        rngval_[i] = 0.;
      } else if (uba[i] > CPX_INFBOUND) {
        // Ineq - no upper bound
        sense_[i] = 'G';
        rhs_[i] = lba[i];
        rngval_[i] = 0.;
      } else { // Inew both upper and lower bounds
        sense_[i] = 'R';
        rhs_[i] = lba[i];
        rngval_[i] = uba[i] - lba[i];
      }
    }

    // Copying objective, constraints, and bounds.
    const Sparsity& A_sp = input(QP_SOLVER_A).sparsity();
    const int* matbeg = getPtr(A_sp.colind());
    const int* matind = getPtr(A_sp.row());
    const double* matval = input(QP_SOLVER_A).ptr();
    const double* obj = input(QP_SOLVER_G).ptr();
    const double* lb = input(QP_SOLVER_LBX).ptr();
    const double* ub = input(QP_SOLVER_UBX).ptr();
    status = CPXcopylp(env_, lp_, n_, nc_, objsen_, obj, rhs_.data(), sense_.data(),
                        matbeg, getPtr(matcnt_), matind, matval, lb, ub, rngval_.data());

    // Preparing coefficient matrix Q
    const Sparsity& H_sp = input(QP_SOLVER_H).sparsity();
    const int* qmatbeg = getPtr(H_sp.colind());
    const int* qmatind = getPtr(H_sp.row());
    const double* qmatval = input(QP_SOLVER_H).ptr();
    status = CPXcopyquad(env_, lp_, qmatbeg, getPtr(qmatcnt_), qmatind, qmatval);

    if (dump_to_file_) {
      const char* fn = string(getOption("dump_filename")).c_str();
      CPXwriteprob(env_, lp_, fn, "LP");
    }

    // Warm-starting if possible
    const double* x0 = input(QP_SOLVER_X0).ptr();
    const double* lam_x0 = input(QP_SOLVER_LAM_X0).ptr();
    if (qp_method_ != 0 && qp_method_ != 4 && is_warm_) {
      // TODO(Joel): Initialize slacks and dual variables of bound constraints
      CPXcopystart(env_, lp_, getPtr(cstat_), getPtr(rstat_), x0, NULL, NULL, lam_x0);
    } else {
      status = CPXcopystart(env_, lp_, NULL, NULL, x0, NULL, NULL, lam_x0);
    }

    // Optimize...
    status = CPXqpopt(env_, lp_);

    if (status) {
      casadi_error("CPLEX: Failed to solve QP...");
    }
    // Retrieving solution
    int solstat;

    std::vector<double> slack;
    slack.resize(nc_);
    status = CPXsolution(env_, lp_, &solstat,
                          output(QP_SOLVER_COST).ptr(),
                          output(QP_SOLVER_X).ptr(),
                          output(QP_SOLVER_LAM_A).ptr(),
                          getPtr(slack),
                          output(QP_SOLVER_LAM_X).ptr());

    if (status) {
      cout << "CPLEX: Failed to get solution.\n";
    }
    // Retrieving the basis
    if (qp_method_ != 0 && qp_method_ != 4) {
      status = CPXgetbase(env_, lp_, getPtr(cstat_), getPtr(rstat_));
    }

    // Flip the sign of the multipliers
    for (vector<double>::iterator it=output(QP_SOLVER_LAM_A).begin();
        it!=output(QP_SOLVER_LAM_A).end(); ++it) *it = -*it;
    for (vector<double>::iterator it=output(QP_SOLVER_LAM_X).begin();
        it!=output(QP_SOLVER_LAM_X).end(); ++it) *it = -*it;

    int solnstat = CPXgetstat(env_, lp_);
    stringstream errormsg; // NOTE: Why not print directly to cout and cerr?
    if (verbose()) {
      if (solnstat == CPX_STAT_OPTIMAL) {
        errormsg << "CPLEX: solution status: Optimal solution found.\n";
      } else if (solnstat == CPX_STAT_UNBOUNDED) {
        errormsg << "CPLEX: solution status: Model is unbounded\n";
      } else if (solnstat == CPX_STAT_INFEASIBLE) {
        errormsg << "CPLEX: solution status: Model is infeasible\n";
      } else if (solnstat == CPX_STAT_INForUNBD) {
        errormsg << "CPLEX: solution status: Model is infeasible or unbounded\n";
      } else if (solnstat == CPX_STAT_OPTIMAL_INFEAS) {
        errormsg << "CPLEX: solution status: Optimal solution "
            "is available but with infeasibilities\n";
      } else if (solnstat == CPX_STAT_NUM_BEST) {
        errormsg << "CPLEX: solution status: Solution available, but not "
            "proved optimal due to numeric difficulties.\n";
      } else if (solnstat == CPX_STAT_FIRSTORDER) {
        errormsg << "CPLEX: solution status: Solution satisfies first-order optimality "
            "conditions, but is not necessarily globally optimal.\n";
      } else {
        errormsg << "CPLEX: solution status: " <<  solnstat << "\n";
      }
      cout << errormsg.str();

      // Printing basis condition number
      //double cn;
      //status = CPXgetdblquality(env_, lp_, &cn, CPX_KAPPA);
      //cout << "CPLEX: Basis condition number: " << cn << endl;
    }
    if (solnstat != CPX_STAT_OPTIMAL) {
      //    throw CasadiException(errormsg.c_str());
    }

    // Next time we warm start
    if (static_cast<bool>(getOption("warm_start"))) {
      is_warm_ = true;
    }

  }

  CplexInterface* CplexInterface::clone() const {
    // Return a deepcopy
    CplexInterface* node = new CplexInterface(st_);
    if (!node->is_init_)
      node->init();
    return node;
  }

  CplexInterface::~CplexInterface() {
    freeCplex();
  }

  void CplexInterface::freeCplex() {
    // Return flag
    int status;

    // Only free if Cplex problem if it has been allocated
    if (lp_!=0) {
      status = CPXfreeprob(env_, &lp_);
      if (status!=0) {
        std::cerr << "CPXfreeprob failed, error code " << status << ".\n";
      }
      lp_ = 0;
    }

    // Closing down license
    if (env_!=0) {
      status = CPXcloseCPLEX(&env_);
      if (status!=0) {
        std::cerr << "CPXcloseCPLEX failed, error code " << status << ".\n";
      }
      env_ = 0;
    }
  }

} // end namespace casadi
