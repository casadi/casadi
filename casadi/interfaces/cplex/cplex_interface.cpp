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
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <vector>

#include "ilcplex/cplex.h"

namespace casadi {

  using namespace std;

  extern "C"
  int CASADI_QPSOL_CPLEX_EXPORT
  casadi_register_qpsol_cplex(Qpsol::Plugin* plugin) {
    plugin->creator = CplexInterface::creator;
    plugin->name = "cplex";
    plugin->doc = CplexInterface::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_QPSOL_CPLEX_EXPORT casadi_load_qpsol_cplex() {
    Qpsol::registerPlugin(casadi_register_qpsol_cplex);
  }

  CplexInterface::CplexInterface(const std::string& name,
                                 const std::map<std::string, Sparsity>& st)
    : Qpsol(name, st) {

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
    Qpsol::init();

    qp_method_     = optionEnumValue("qp_method");
    dump_to_file_  = option("dump_to_file");
    tol_           = option("tol");
    //  dump_filename_ = option("dump_filename");

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
      userOut() << "CPLEX: Problem with setting parameter... ERROR: " << status << std::endl;
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
    status = CPXsetintparam(env_, CPX_PARAM_DEPIND, optionEnumValue("dep_check"));
    // Setting barrier iteration limit
    status = CPXsetintparam(env_, CPX_PARAM_BARITLIM, option("barrier_maxiter"));
    // Setting simplex iteration limit
    status = CPXsetintparam(env_, CPX_PARAM_ITLIM, option("simplex_maxiter"));
    if (qp_method_ == 7) {
      // Setting crossover algorithm
      status = CPXsetintparam(env_, CPX_PARAM_BARCROSSALG, 1);
    }
    if (!static_cast<bool>(option("convex"))) {
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
    const Sparsity& A_sp = input(QPSOL_A).sparsity();
    matcnt_.resize(A_sp.size2());
    transform(A_sp.colind()+1, A_sp.colind() + A_sp.size2()+1, A_sp.colind(), matcnt_.begin(),
              minus<int>());

    // Matrix H, count the number of elements per column
    const Sparsity& H_sp = input(QPSOL_H).sparsity();
    qmatcnt_.resize(H_sp.size2());
    transform(H_sp.colind()+1, H_sp.colind() + H_sp.size2()+1, H_sp.colind(), qmatcnt_.begin(),
              minus<int>());

    casadi_assert(lp_==0);
    lp_ = CPXcreateprob(env_, &status, "QP from CasADi");
  }

  void CplexInterface::evalD(void* mem, const double** arg, double** res, int* iw, double* w) {
    // Number of inputs and outputs
    int num_in = n_in();
    int num_out = n_out();

    // Pass the inputs to the function
    for (int i=0; i<num_in; ++i) {
      if (arg[i] != 0) {
        setInputNZ(arg[i], i);
      } else {
        setInput(0., i);
      }
    }

    if (inputs_check_) checkInputs();

    int status;

    // We change method in crossover
    if (is_warm_ && qp_method_ == 7) {
      status = CPXsetintparam(env_, CPX_PARAM_QPMETHOD, 1);
    }

    // Looping over constraints
    const vector<double>& lba = input(QPSOL_LBA).data();
    const vector<double>& uba = input(QPSOL_UBA).data();

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
    const Sparsity& A_sp = input(QPSOL_A).sparsity();
    const int* matbeg = A_sp.colind();
    const int* matind = A_sp.row();
    const double* matval = input(QPSOL_A).ptr();
    const double* obj = input(QPSOL_G).ptr();
    const double* lb = input(QPSOL_LBX).ptr();
    const double* ub = input(QPSOL_UBX).ptr();
    status = CPXcopylp(env_, lp_, n_, nc_, objsen_, obj, rhs_.data(), sense_.data(),
                        matbeg, getPtr(matcnt_), matind, matval, lb, ub, rngval_.data());

    // Preparing coefficient matrix Q
    const Sparsity& H_sp = input(QPSOL_H).sparsity();
    const int* qmatbeg = H_sp.colind();
    const int* qmatind = H_sp.row();
    const double* qmatval = input(QPSOL_H).ptr();
    status = CPXcopyquad(env_, lp_, qmatbeg, getPtr(qmatcnt_), qmatind, qmatval);

    if (dump_to_file_) {
      const char* fn = string(option("dump_filename")).c_str();
      CPXwriteprob(env_, lp_, fn, "LP");
    }

    // Warm-starting if possible
    const double* x0 = input(QPSOL_X0).ptr();
    const double* lam_x0 = input(QPSOL_LAM_X0).ptr();
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
                          output(QPSOL_COST).ptr(),
                          output(QPSOL_X).ptr(),
                          output(QPSOL_LAM_A).ptr(),
                          getPtr(slack),
                          output(QPSOL_LAM_X).ptr());

    if (status) {
      userOut() << "CPLEX: Failed to get solution.\n";
    }
    // Retrieving the basis
    if (qp_method_ != 0 && qp_method_ != 4) {
      status = CPXgetbase(env_, lp_, getPtr(cstat_), getPtr(rstat_));
    }

    // Flip the sign of the multipliers
    for (auto it=output(QPSOL_LAM_A)->begin();
         it!=output(QPSOL_LAM_A)->end(); ++it) *it = -*it;
    for (auto it=output(QPSOL_LAM_X)->begin();
         it!=output(QPSOL_LAM_X)->end(); ++it) *it = -*it;

    int solnstat = CPXgetstat(env_, lp_);
    stringstream errormsg;
    // NOTE: Why not print directly to userOut() and userOut<true, PL_WARN>()?
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
      userOut() << errormsg.str();

      // Printing basis condition number
      //double cn;
      //status = CPXgetdblquality(env_, lp_, &cn, CPX_KAPPA);
      //userOut() << "CPLEX: Basis condition number: " << cn << endl;
    }
    if (solnstat != CPX_STAT_OPTIMAL) {
      //    throw CasadiException(errormsg.c_str());
    }

    // Next time we warm start
    if (static_cast<bool>(option("warm_start"))) {
      is_warm_ = true;
    }

    // Get the outputs
    for (int i=0; i<num_out; ++i) {
      if (res[i] != 0) getOutputNZ(res[i], i);
    }
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
        userOut<true, PL_WARN>() << "CPXfreeprob failed, error code " << status << ".\n";
      }
      lp_ = 0;
    }

    // Closing down license
    if (env_!=0) {
      status = CPXcloseCPLEX(&env_);
      if (status!=0) {
        userOut<true, PL_WARN>() << "CPXcloseCPLEX failed, error code " << status << ".\n";
      }
      env_ = 0;
    }
  }

} // end namespace casadi
