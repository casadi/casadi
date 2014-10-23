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


#include "worhp_interface.hpp"
#include "casadi/core/std_vector_tools.hpp"
#include "casadi/core/matrix/matrix_tools.hpp"
#include "casadi/core/mx/mx_tools.hpp"
#include "casadi/core/matrix/sparsity_tools.hpp"
#include "casadi/core/function/mx_function.hpp"
#include <ctime>

using namespace std;

namespace casadi {

  extern "C"
  int CASADI_NLPSOLVER_WORHP_EXPORT
  casadi_register_nlpsolver_worhp(NlpSolverInternal::Plugin* plugin) {
    plugin->creator = WorhpInterface::creator;
    plugin->name = "worhp";
    plugin->doc = WorhpInterface::meta_doc.c_str();
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOLVER_WORHP_EXPORT casadi_load_nlpsolver_worhp() {
    NlpSolverInternal::registerPlugin(casadi_register_nlpsolver_worhp);
  }

  WorhpInterface::WorhpInterface(const Function& nlp) : NlpSolverInternal(nlp) {

    // Monitors
    addOption("monitor",            OT_STRINGVECTOR,  GenericType(),  "Monitor functions",
              "eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h", true);
    addOption("print_time",         OT_BOOLEAN,       true,
              "Print information about execution time");

    int status;
    InitParams(&status, &worhp_p_);

    std::stringstream ss;
    ss << "Armijo recovery strategies. Vector of size " << NAres;

    std::vector<int> ares(NAres);
    std::copy(worhp_p_.Ares, worhp_p_.Ares+NAres, ares.begin());
    addOption("Ares", OT_INTEGERVECTOR, ares, ss.str());

    for (int i=0;i<WorhpGetParamCount();++i) {
      WorhpType type = WorhpGetParamType(i+1);
      const char* name = WorhpGetParamName(i+1);
      if (strcmp(name, "Ares")==0) continue;
      switch (type) {
        case WORHP_BOOL_T:
          bool default_bool;
          WorhpGetBoolParam(&worhp_p_, name, &default_bool);
          addOption(WorhpGetParamName(i+1), OT_BOOLEAN, default_bool,
                    WorhpGetParamDescription(i+1));
          break;
        case WORHP_DOUBLE_T:
          double default_double;
          WorhpGetDoubleParam(&worhp_p_, name, &default_double);
          addOption(WorhpGetParamName(i+1), OT_REAL, default_double, WorhpGetParamDescription(i+1));
          break;
        case WORHP_INT_T:
          int default_int;
          WorhpGetIntParam(&worhp_p_, name, &default_int);
          addOption(WorhpGetParamName(i+1), OT_INTEGER, default_int, WorhpGetParamDescription(i+1));
          break;
        default:
          break;// do nothing
      }
    }

    addOption("qp_ipBarrier", OT_REAL, worhp_p_.qp.ipBarrier, "IP barrier parameter.");
    addOption("qp_ipComTol", OT_REAL, worhp_p_.qp.ipComTol, "IP complementarity tolerance.");
    addOption("qp_ipFracBound", OT_REAL, worhp_p_.qp.ipFracBound,
              "IP fraction-to-the-boundary parameter.");
    addOption("qp_ipLsMethod", OT_STRING, GenericType(),
      "Select the direct linear solver used by the IP method.",
      "LAPACK::0|MA57: only available if provided by the user:1|"
      "SuperLU::2|PARDISO: only available if provided by the user, "
      "subject to license availability:3|"
      "MUMPS: currently Linux platforms only:5|"
      "WSMP: subject to license availability:6|"
      "MA86: experimental, only available if provided by the user:7|"
      "MA97:experimental, only available if provided by the user:8");
    setOptionByEnumValue("qp_ipLsMethod", worhp_p_.qp.ipLsMethod);
    addOption("qp_ipMinAlpha", OT_REAL, worhp_p_.qp.ipMinAlpha,
              "IP line search minimum step size.");
    addOption("qp_ipTryRelax", OT_BOOLEAN, worhp_p_.qp.ipTryRelax,
      "Enable relaxation strategy when encountering an error.");
    addOption("qp_ipRelaxDiv", OT_REAL, worhp_p_.qp.ipRelaxDiv,
      "The relaxation term is divided by this value if successful.");
    addOption("qp_ipRelaxMult", OT_REAL, worhp_p_.qp.ipRelaxMult,
      "The relaxation term is multiplied by this value if unsuccessful.");
    addOption("qp_ipRelaxMax", OT_REAL, worhp_p_.qp.ipRelaxMax, "Maximum relaxation value.");
    addOption("qp_ipRelaxMin", OT_REAL, worhp_p_.qp.ipRelaxMin, "Mimimum relaxation value.");
    addOption("qp_ipResTol", OT_REAL, worhp_p_.qp.ipResTol, "IP residuals tolerance.");
    addOption("qp_lsItMaxIter", OT_INTEGER, worhp_p_.qp.lsItMaxIter,
      "Maximum number of iterations of the iterative linear solvers.");
    addOption("qp_lsItMethod", OT_STRING, GenericType(),
      "Select the iterative linear solver.",
      "none:Deactivate; use a direct linear solver.:0|CGNR::1|CGNE::2|CGS::3|BiCGSTAB::4");
    setOptionByEnumValue("qp_lsItMethod", worhp_p_.qp.lsItMethod);
    addOption("qp_lsItPrecondMethod", OT_STRING, GenericType(),
      "Select preconditioner for the iterative linear solver.",
      "none:No preconditioner.:0|"
      "static:Static preconditioner (KKT-matrix with constant lower-right block).:1|"
      "full:Full KKT-matrix.:2");
    setOptionByEnumValue("qp_lsItPrecondMethod", worhp_p_.qp.lsItPrecondMethod);
    addOption("qp_lsRefineMaxIter", OT_INTEGER, worhp_p_.qp.lsRefineMaxIter,
      "Maximum number of iterative refinement steps of the direct linear solvers.");
    addOption("qp_lsScale", OT_BOOLEAN, worhp_p_.qp.lsScale,
              "Enables scaling on linear solver level.");
    addOption("qp_lsTrySimple", OT_BOOLEAN, worhp_p_.qp.lsTrySimple,
      "Some matrices can be solved without calling a linear equation solver."
      "Currently only diagonal matrices are supported."
      "Non-diagonal matrices will besolved with the chosen linear equation solver.");
    addOption("qp_lsTol", OT_REAL, worhp_p_.qp.lsTol, "Tolerance for the linear solver.");
    addOption("qp_maxIter", OT_INTEGER, worhp_p_.qp.maxIter,
      "Imposes an upper limit on the number of minor solver iterations, "
      " i.e. for the quadratic subproblem solver."
      "If the limit is reached before convergence, "
      "WORHP will activate QP recovery strategies to prevent a solver breakdown.");
    addOption("qp_method", OT_STRING, GenericType(),
      "Select the solution method used by the QP solver.",
      "ip:Interior-Point method.:1|nsn:Nonsmooth-Newton method.:2|"
      "automatic: Prefer IP and fall back to NSN on error.:12");
    setOptionByEnumValue("qp_method", worhp_p_.qp.method);
    addOption("qp_nsnBeta", OT_REAL, worhp_p_.qp.nsnBeta, "NSN stepsize decrease factor.");
    addOption("qp_nsnGradStep", OT_BOOLEAN, worhp_p_.qp.nsnGradStep,
      "Enable gradient steps in the NSN method.");
    addOption("qp_nsnKKT", OT_REAL, worhp_p_.qp.nsnKKT, "NSN KKT tolerance.");
    addOption("qp_nsnLsMethod", OT_STRING, GenericType(),
      "Select the direct linear solver used by the NSN method.",
      "SuperLU::2|MA48: only available if provided by the user:4");
    setOptionByEnumValue("qp_nsnLsMethod", worhp_p_.qp.nsnLsMethod);
    addOption("qp_nsnMinAlpha", OT_REAL, worhp_p_.qp.nsnMinAlpha,
              "NSN line search minimum step size.");
    addOption("qp_nsnSigma", OT_REAL, worhp_p_.qp.nsnSigma, "NSN line search slope parameter.");
    addOption("qp_printLevel", OT_STRING, GenericType(),
      "Controls the amount of QP solver output.",
      "none:No output.:0|warn:Print warnings and errors.:1|iterations:Print iterations.:2");
    setOptionByEnumValue("qp_printLevel", worhp_p_.qp.printLevel);
    addOption("qp_scaleIntern", OT_BOOLEAN, worhp_p_.qp.scaleIntern, "Enable scaling on QP level.");
    addOption("qp_strict", OT_BOOLEAN, worhp_p_.qp.strict,
              "Use strict termination criteria in IP method.");

    worhp_o_.initialised = false;
    worhp_w_.initialised = false;
    worhp_p_.initialised = false;
    worhp_c_.initialised = false;

    // WORKAROUND: Bug in scaling, set to false by default // FIXME
    setOption("ScaledObj", false);

    // WORKAROUND: Why is this needed? // FIXME
    setOption("ScaleConIter", true);
  }

  WorhpInterface::~WorhpInterface() {
    if (worhp_p_.initialised || worhp_o_.initialised ||
        worhp_w_.initialised || worhp_c_.initialised)
      WorhpFree(&worhp_o_, &worhp_w_, &worhp_p_, &worhp_c_);
  }

  WorhpInterface* WorhpInterface::clone() const {
    // Use default copy routine
    WorhpInterface* node = new WorhpInterface(*this);

    // Mark Worhp datastructures not initialized to avoid double freeing
    node->worhp_o_.initialised = false;
    node->worhp_w_.initialised = false;
    node->worhp_p_.initialised = false;
    node->worhp_c_.initialised = false;

    return node;
  }

  void WorhpInterface::init() {

    // Call the init method of the base class
    NlpSolverInternal::init();

    if (hasSetOption("Ares")) {
      std::vector<int> ares = getOption("Ares");
      std::copy(ares.begin(), ares.begin()+NAres, worhp_p_.Ares);
    }

    // Read options
    passOptions();

    // Exact Hessian?
    exact_hessian_ = getOption("UserHM");

    // Get/generate required functions
    gradF();
    jacG();
    if (exact_hessian_) { // does not appear to work
      hessLag();
    }

    // Update status?
    status_[TerminateSuccess]="TerminateSuccess";
    status_[OptimalSolution]="OptimalSolution";
    status_[SearchDirectionZero]="SearchDirectionZero";
    status_[SearchDirectionSmall]="SearchDirectionSmall";
    status_[StationaryPointFound]="StationaryPointFound";
    status_[AcceptableSolution]="AcceptableSolution";
    status_[AcceptablePrevious]="AcceptablePrevious";
    status_[FritzJohn]="FritzJohn";
    status_[NotDiffable]="NotDiffable";
    status_[Unbounded]="Unbounded";
    status_[FeasibleSolution]="FeasibleSolution";
    status_[LowPassFilterOptimal]="LowPassFilterOptimal";
    status_[LowPassFilterAcceptable]="LowPassFilterAcceptable";
    status_[TerminateError]="TerminateError";
    status_[InitError]="InitError";
    status_[DataError]="DataError";
    status_[MaxCalls]="MaxCalls";
    status_[MaxIter]="MaxIter";
    status_[MinimumStepsize]="MinimumStepsize";
    status_[QPerror]="QPerror";
    status_[ProblemInfeasible]="ProblemInfeasible";
    status_[GroupsComposition]="GroupsComposition";
    status_[TooBig]="TooBig";
    status_[Timeout]="Timeout";
    status_[FDError]="FDError";
    status_[LocalInfeas]="LocalInfeas";
    status_[LicenseError]="LicenseError. Please set the WORHP_LICENSE_FILE environmental "
        "variable with the full path to the license file";
    status_[TerminatedByUser]="TerminatedByUser";
    status_[FunctionErrorF]="FunctionErrorF";
    status_[FunctionErrorG]="FunctionErrorG";
    status_[FunctionErrorDF]="FunctionErrorDF";
    status_[FunctionErrorDG]="FunctionErrorDG";
    status_[FunctionErrorHM]="FunctionErrorHM";
  }

  void WorhpInterface::reset() {

    // Number of (free) variables
    worhp_o_.n = nx_;

    // Number of constraints
    worhp_o_.m = ng_;

    // Free existing Worhp memory (except parameters)
    bool p_init_backup = worhp_p_.initialised;
    worhp_p_.initialised = false; // Avoid freeing the memory for parameters
    if (worhp_o_.initialised || worhp_w_.initialised || worhp_c_.initialised) {
      WorhpFree(&worhp_o_, &worhp_w_, &worhp_p_, &worhp_c_);
    }
    worhp_p_.initialised = p_init_backup;

    /// Control data structure needs to be reset every time
    worhp_c_.initialised = false;
    worhp_w_.initialised = false;
    worhp_o_.initialised = false;

    // Worhp uses the CS format internally, hence it is the preferred sparse matrix format.
    worhp_w_.DF.nnz = nx_;
    if (worhp_o_.m>0) {
      worhp_w_.DG.nnz = jacG_.output().size();  // Jacobian of G
    } else {
      worhp_w_.DG.nnz = 0;
    }

    if (exact_hessian_ /*worhp_w_.HM.NeedStructure*/) { // not initialized

      // Get the sparsity pattern of the Hessian
      const Sparsity& spHessLag = this->spHessLag();
      const vector<int>& colind = spHessLag.colind();
      const vector<int>& row = spHessLag.row();

      // Get number of nonzeros in the lower triangular part of the Hessian including full diagonal
      worhp_w_.HM.nnz = nx_; // diagonal entries
      for (int c=0; c<nx_; ++c) {
        for (int el=colind[c]; el<colind[c+1] && row[el]<c; ++el) {
          worhp_w_.HM.nnz++; // strictly lower triangular part
        }
      }
    } else {
      worhp_w_.HM.nnz = 0;
    }

    /* Data structure initialisation. */
    WorhpInit(&worhp_o_, &worhp_w_, &worhp_p_, &worhp_c_);
    if (worhp_c_.status != FirstCall) {
      casadi_error("Main: Initialisation failed. Status: " << formatStatus(worhp_c_.status));
    }

    if (worhp_w_.DF.NeedStructure) {
      for (int i=0; i<nx_; ++i) {
        worhp_w_.DF.row[i] = i + 1; // Index-1 based
      }
    }

    if (worhp_o_.m>0 && worhp_w_.DG.NeedStructure) {
      // Get sparsity pattern. Note WORHP is column major
      const DMatrix & J = jacG_.output(JACG_JAC);

      int nz=0;
      const vector<int>& colind = J.colind();
      const vector<int>& row = J.row();
      for (int c=0; c<nx_; ++c) {
        for (int el=colind[c]; el<colind[c+1]; ++el) {
          int r = row[el];
          worhp_w_.DG.col[nz] = c + 1; // Index-1 based
          worhp_w_.DG.row[nz] = r + 1;
          nz++;
        }
      }
    }

    if (worhp_w_.HM.NeedStructure) {
      // Get the sparsity pattern of the Hessian
      const Sparsity& spHessLag = this->spHessLag();
      const vector<int>& colind = spHessLag.colind();
      const vector<int>& row = spHessLag.row();

      int nz=0;

      // Upper triangular part of the Hessian (note CCS -> CRS format change)
      for (int c=0; c<nx_; ++c) {
        for (int el=colind[c]; el<colind[c+1]; ++el) {
          if (row[el]>c) {
            worhp_w_.HM.row[nz] = row[el] + 1;
            worhp_w_.HM.col[nz] = c + 1;
            nz++;
          }
        }
      }

      // Diagonal always included
      for (int r=0; r<nx_; ++r) {
        worhp_w_.HM.row[nz] = r + 1;
        worhp_w_.HM.col[nz] = r + 1;
        nz++;
      }
    }
  }

  void WorhpInterface::setQPOptions() {
    setOption("UserHM", true);
  }

  void WorhpInterface::passOptions() {

     for (int i=0;i<WorhpGetParamCount();++i) {
      WorhpType type = WorhpGetParamType(i+1);
      const char* name = WorhpGetParamName(i+1);
      if (strcmp(name, "Ares")==0) continue;

      switch (type) {
        case WORHP_BOOL_T:
          if (hasSetOption(name)) WorhpSetBoolParam(&worhp_p_, name, getOption(name));
          break;
        case WORHP_DOUBLE_T:
          if (hasSetOption(name)) WorhpSetDoubleParam(&worhp_p_, name, getOption(name));
          break;
        case WORHP_INT_T:
          if (hasSetOption(name)) WorhpSetIntParam(&worhp_p_, name, getOption(name));
          break;
        default:
          break;// do nothing
      }
    }

    if (hasSetOption("qp_ipBarrier")) worhp_p_.qp.ipBarrier = getOption("qp_ipBarrier");
    if (hasSetOption("qp_ipComTol")) worhp_p_.qp.ipComTol = getOption("qp_ipComTol");
    if (hasSetOption("qp_ipFracBound")) worhp_p_.qp.ipFracBound = getOption("qp_ipFracBound");
    if (hasSetOption("qp_ipLsMethod")) worhp_p_.qp.ipLsMethod = getOptionEnumValue("qp_ipLsMethod");
    if (hasSetOption("qp_ipMinAlpha")) worhp_p_.qp.ipMinAlpha = getOption("qp_ipMinAlpha");
    if (hasSetOption("qp_ipTryRelax")) worhp_p_.qp.ipTryRelax = getOption("qp_ipTryRelax");
    if (hasSetOption("qp_ipRelaxDiv")) worhp_p_.qp.ipRelaxDiv = getOption("qp_ipRelaxDiv");
    if (hasSetOption("qp_ipRelaxMult")) worhp_p_.qp.ipRelaxMult = getOption("qp_ipRelaxMult");
    if (hasSetOption("qp_ipRelaxMax")) worhp_p_.qp.ipRelaxMax = getOption("qp_ipRelaxMax");
    if (hasSetOption("qp_ipRelaxMin")) worhp_p_.qp.ipRelaxMin = getOption("qp_ipRelaxMin");
    if (hasSetOption("qp_ipResTol")) worhp_p_.qp.ipResTol = getOption("qp_ipResTol");
    if (hasSetOption("qp_lsItMaxIter")) worhp_p_.qp.lsItMaxIter = getOption("qp_lsItMaxIter");
    if (hasSetOption("qp_lsItMethod")) worhp_p_.qp.lsItMethod = getOptionEnumValue("qp_lsItMethod");
    if (hasSetOption("qp_lsItPrecondMethod"))
        worhp_p_.qp.lsItPrecondMethod = getOptionEnumValue("qp_lsItPrecondMethod");
    if (hasSetOption("qp_lsRefineMaxIter"))
        worhp_p_.qp.lsRefineMaxIter = getOption("qp_lsRefineMaxIter");
    if (hasSetOption("qp_lsScale")) worhp_p_.qp.lsScale = getOption("qp_lsScale");
    if (hasSetOption("qp_lsTrySimple")) worhp_p_.qp.lsTrySimple = getOption("qp_lsTrySimple");
    if (hasSetOption("qp_lsTol")) worhp_p_.qp.lsTol = getOption("qp_lsTol");
    if (hasSetOption("qp_maxIter")) worhp_p_.qp.maxIter = getOption("qp_maxIter");
    if (hasSetOption("qp_method")) worhp_p_.qp.method = getOptionEnumValue("qp_method");
    if (hasSetOption("qp_nsnBeta")) worhp_p_.qp.nsnBeta = getOption("qp_nsnBeta");
    if (hasSetOption("qp_nsnGradStep")) worhp_p_.qp.nsnGradStep = getOption("qp_nsnGradStep");
    if (hasSetOption("qp_nsnKKT")) worhp_p_.qp.nsnKKT = getOption("qp_nsnKKT");
    if (hasSetOption("qp_nsnLsMethod"))
        worhp_p_.qp.nsnLsMethod = getOptionEnumValue("qp_nsnLsMethod");
    if (hasSetOption("qp_nsnMinAlpha")) worhp_p_.qp.nsnMinAlpha = getOption("qp_nsnMinAlpha");
    if (hasSetOption("qp_nsnSigma")) worhp_p_.qp.nsnSigma = getOption("qp_nsnSigma");
    if (hasSetOption("qp_printLevel")) worhp_p_.qp.printLevel = getOptionEnumValue("qp_printLevel");
    if (hasSetOption("qp_scaleIntern")) worhp_p_.qp.scaleIntern = getOption("qp_scaleIntern");
    if (hasSetOption("qp_strict")) worhp_p_.qp.strict = getOption("qp_strict");

    // Mark the parameters as set
    worhp_p_.initialised = true;
  }

  std::string WorhpInterface::formatStatus(int status) const {
    if (status_.find(status)==status_.end()) {
      std::stringstream ss;
      ss << "Unknown status: " << status;
      return ss.str();
    } else {
      return (*status_.find(status)).second;
    }
  }

  void WorhpInterface::evaluate() {
    log("WorhpInterface::evaluate");

    if (gather_stats_) {
      Dictionary iterations;
      iterations["iter_sqp"] = std::vector<int>();
      iterations["inf_pr"] = std::vector<double>();
      iterations["inf_du"] = std::vector<double>();
      iterations["obj"] = std::vector<double>();
      iterations["alpha_pr"] = std::vector<double>();
      stats_["iterations"] = iterations;
    }

    // Prepare the solver
    reset();

    if (inputs_check_) checkInputs();
    checkInitialBounds();

    // Reset the counters
    t_eval_f_ = t_eval_grad_f_ = t_eval_g_ = t_eval_jac_g_ = t_eval_h_ = t_callback_fun_ =
      t_callback_prepare_ = t_mainloop_ = 0;

    n_eval_f_ = n_eval_grad_f_ = n_eval_g_ = n_eval_jac_g_ = n_eval_h_ = 0;

    // Get inputs
    log("WorhpInterface::evaluate: Reading user inputs");
    const DMatrix& x0 = input(NLP_SOLVER_X0);
    const DMatrix& lbx = input(NLP_SOLVER_LBX);
    const DMatrix& ubx = input(NLP_SOLVER_UBX);
    const DMatrix& lam_x0 = input(NLP_SOLVER_LAM_X0);
    const DMatrix& lbg = input(NLP_SOLVER_LBG);
    const DMatrix& ubg = input(NLP_SOLVER_UBG);
    const DMatrix& lam_g0 = input(NLP_SOLVER_LAM_G0);

    double inf = numeric_limits<double>::infinity();

    for (int i=0;i<nx_;++i) {
      casadi_assert_message(lbx.at(i)!=ubx.at(i),
      "WorhpInterface::evaluate: Worhp cannot handle the case when LBX == UBX."
      "You have that case at non-zero " << i << " , which has value " << ubx.at(i) << "."
      "Reformulate your problem by using a parameter for the corresponding variable.");
    }

    for (int i=0;i<lbg.size();++i) {
      casadi_assert_message(!(lbg.at(i)==-inf && ubg.at(i) == inf),
      "WorhpInterface::evaluate: Worhp cannot handle the case when both LBG and UBG are infinite."
      "You have that case at non-zero " << i << "."
      "Reformulate your problem eliminating the corresponding constraint.");
    }

    // Pass inputs to WORHP data structures
    x0.getArray(worhp_o_.X, worhp_o_.n);
    lbx.getArray(worhp_o_.XL, worhp_o_.n);
    ubx.getArray(worhp_o_.XU, worhp_o_.n);
    lam_x0.getArray(worhp_o_.Lambda, worhp_o_.n);
    if (worhp_o_.m>0) {
      lam_g0.getArray(worhp_o_.Mu, worhp_o_.m);
      lbg.getArray(worhp_o_.GL, worhp_o_.m);
      ubg.getArray(worhp_o_.GU, worhp_o_.m);
    }

    // Replace infinite bounds with worhp_p_.Infty
    for (int i=0; i<nx_; ++i) if (worhp_o_.XL[i]==-inf) worhp_o_.XL[i] = -worhp_p_.Infty;
    for (int i=0; i<nx_; ++i) if (worhp_o_.XU[i]== inf) worhp_o_.XU[i] =  worhp_p_.Infty;
    for (int i=0; i<ng_; ++i) if (worhp_o_.GL[i]==-inf) worhp_o_.GL[i] = -worhp_p_.Infty;
    for (int i=0; i<ng_; ++i) if (worhp_o_.GU[i]== inf) worhp_o_.GU[i] =  worhp_p_.Infty;

    log("WorhpInterface::starting iteration");

    double time1 = clock();

    // Reverse Communication loop
    while (worhp_c_.status < TerminateSuccess &&  worhp_c_.status > TerminateError) {
      if (GetUserAction(&worhp_c_, callWorhp)) {
        Worhp(&worhp_o_, &worhp_w_, &worhp_p_, &worhp_c_);
      }


      if (GetUserAction(&worhp_c_, iterOutput)) {

        if (!worhp_w_.FirstIteration) {
          if (gather_stats_) {
            Dictionary & iterations = stats_["iterations"];
            static_cast<std::vector<int> &>(iterations["iter_sqp"]).push_back(worhp_w_.MinorIter);
            static_cast<std::vector<double> &>(iterations["inf_pr"]).push_back(worhp_w_.NormMax_CV);
            static_cast<std::vector<double> &>(iterations["inf_du"]).push_back(worhp_w_.ScaledKKT);
            static_cast<std::vector<double> &>(iterations["obj"]).push_back(worhp_o_.F);
            static_cast<std::vector<double> &>(
              iterations["alpha_pr"]).push_back(worhp_w_.ArmijoAlpha);
          }

          if (!callback_.isNull()) {
            double time1 = clock();
            // Copy outputs
            if (!output(NLP_SOLVER_X).isEmpty())
              output(NLP_SOLVER_X).setArray(worhp_o_.X, worhp_o_.n);
            if (!output(NLP_SOLVER_F).isEmpty())
              output(NLP_SOLVER_F).set(worhp_o_.F);
            if (!output(NLP_SOLVER_G).isEmpty())
              output(NLP_SOLVER_G).setArray(worhp_o_.G, worhp_o_.m);
            if (!output(NLP_SOLVER_LAM_X).isEmpty())
              output(NLP_SOLVER_LAM_X).setArray(worhp_o_.Lambda, worhp_o_.n);
            if (!output(NLP_SOLVER_LAM_G).isEmpty())
              output(NLP_SOLVER_LAM_G).setArray(worhp_o_.Mu, worhp_o_.m);

            Dictionary iteration;
            iteration["iter"] = worhp_w_.MajorIter;
            iteration["iter_sqp"] = worhp_w_.MinorIter;
            iteration["inf_pr"] = worhp_w_.NormMax_CV;
            iteration["inf_du"] = worhp_w_.ScaledKKT;
            iteration["obj"] = worhp_o_.F;
            iteration["alpha_pr"] = worhp_w_.ArmijoAlpha;
            stats_["iteration"] = iteration;

            double time2 = clock();
            t_callback_prepare_ += (time2-time1)/CLOCKS_PER_SEC;
            time1 = clock();
            int ret = callback_(ref_, user_data_);
            time2 = clock();
            t_callback_fun_ += (time2-time1)/CLOCKS_PER_SEC;

            if (ret) worhp_c_.status = TerminatedByUser;

          }
        }


        IterationOutput(&worhp_o_, &worhp_w_, &worhp_p_, &worhp_c_);
        DoneUserAction(&worhp_c_, iterOutput);
      }

      if (GetUserAction(&worhp_c_, evalF)) {
        eval_f(worhp_o_.X, worhp_w_.ScaleObj, worhp_o_.F);
        DoneUserAction(&worhp_c_, evalF);
      }

      if (GetUserAction(&worhp_c_, evalG)) {
        eval_g(worhp_o_.X, worhp_o_.G);
        DoneUserAction(&worhp_c_, evalG);
      }

      if (GetUserAction(&worhp_c_, evalDF)) {
        eval_grad_f(worhp_o_.X, worhp_w_.ScaleObj, worhp_w_.DF.val);
        DoneUserAction(&worhp_c_, evalDF);
      }

      if (GetUserAction(&worhp_c_, evalDG)) {
        eval_jac_g(worhp_o_.X, worhp_w_.DG.val);
        DoneUserAction(&worhp_c_, evalDG);
      }

      if (GetUserAction(&worhp_c_, evalHM)) {
        eval_h(worhp_o_.X, worhp_w_.ScaleObj, worhp_o_.Mu, worhp_w_.HM.val);
        DoneUserAction(&worhp_c_, evalHM);
      }

      if (GetUserAction(&worhp_c_, fidif)) {
        WorhpFidif(&worhp_o_, &worhp_w_, &worhp_p_, &worhp_c_);
      }
    }

    double time2 = clock();
    t_mainloop_ += (time2-time1)/CLOCKS_PER_SEC;

    // Copy outputs
    output(NLP_SOLVER_X).setArray(worhp_o_.X, worhp_o_.n, DENSE);
    output(NLP_SOLVER_F).set(worhp_o_.F);
    output(NLP_SOLVER_G).setArray(worhp_o_.G, worhp_o_.m, DENSE);
    output(NLP_SOLVER_LAM_X).setArray(worhp_o_.Lambda, worhp_o_.n);
    output(NLP_SOLVER_LAM_G).setArray(worhp_o_.Mu, worhp_o_.m, DENSE);

    StatusMsg(&worhp_o_, &worhp_w_, &worhp_p_, &worhp_c_);

    if (hasOption("print_time") && static_cast<bool>(getOption("print_time"))) {
      // Write timings
      cout << "time spent in eval_f: " << t_eval_f_ << " s.";
      if (n_eval_f_>0)
        cout << " (" << n_eval_f_ << " calls, " << (t_eval_f_/n_eval_f_)*1000 << " ms. average)";
      cout << endl;
      cout << "time spent in eval_grad_f: " << t_eval_grad_f_ << " s.";
      if (n_eval_grad_f_>0)
        cout << " (" << n_eval_grad_f_ << " calls, "
             << (t_eval_grad_f_/n_eval_grad_f_)*1000 << " ms. average)";
      cout << endl;
      cout << "time spent in eval_g: " << t_eval_g_ << " s.";
      if (n_eval_g_>0)
        cout << " (" << n_eval_g_ << " calls, " << (t_eval_g_/n_eval_g_)*1000 << " ms. average)";
      cout << endl;
      cout << "time spent in eval_jac_g: " << t_eval_jac_g_ << " s.";
      if (n_eval_jac_g_>0)
        cout << " (" << n_eval_jac_g_ << " calls, "
             << (t_eval_jac_g_/n_eval_jac_g_)*1000 << " ms. average)";
      cout << endl;
      cout << "time spent in eval_h: " << t_eval_h_ << " s.";
      if (n_eval_h_>1)
        cout << " (" << n_eval_h_ << " calls, " << (t_eval_h_/n_eval_h_)*1000 << " ms. average)";
      cout << endl;
      cout << "time spent in main loop: " << t_mainloop_ << " s." << endl;
      cout << "time spent in callback function: " << t_callback_fun_ << " s." << endl;
      cout << "time spent in callback preparation: " << t_callback_prepare_ << " s." << endl;
    }

    stats_["t_eval_f"] = t_eval_f_;
    stats_["t_eval_grad_f"] = t_eval_grad_f_;
    stats_["t_eval_g"] = t_eval_g_;
    stats_["t_eval_jac_g"] = t_eval_jac_g_;
    stats_["t_eval_h"] = t_eval_h_;
    stats_["t_mainloop"] = t_mainloop_;
    stats_["t_callback_fun"] = t_callback_fun_;
    stats_["t_callback_prepare"] = t_callback_prepare_;
    stats_["n_eval_f"] = n_eval_f_;
    stats_["n_eval_grad_f"] = n_eval_grad_f_;
    stats_["n_eval_g"] = n_eval_g_;
    stats_["n_eval_jac_g"] = n_eval_jac_g_;
    stats_["n_eval_h"] = n_eval_h_;
    stats_["iter_count"] = worhp_w_.MajorIter;

    stats_["return_code"] = worhp_c_.status;
    stats_["return_status"] = flagmap[worhp_c_.status];

  }

  bool WorhpInterface::eval_h(const double* x, double obj_factor,
                             const double* lambda, double* values) {
    try {
      log("eval_h started");
      double time1 = clock();

      // Make sure generated
      casadi_assert_warning(!hessLag_.isNull(), "Hessian function not pregenerated");

      // Get Hessian
      Function& hessLag = this->hessLag();

      // Pass input
      hessLag.setInput(x, HESSLAG_X);
      hessLag.setInput(input(NLP_SOLVER_P), HESSLAG_P);
      hessLag.setInput(obj_factor, HESSLAG_LAM_F);
      hessLag.setInput(lambda, HESSLAG_LAM_G);

      // Evaluate
      hessLag.evaluate();

      // Get results
      const DMatrix& H = hessLag.output();
      const vector<int>& colind = H.colind();
      const vector<int>& row = H.row();
      const vector<double>& data = H.data();

      // The Hessian values are divided into strictly upper (in WORHP lower) triangular and diagonal
      double* values_upper = values;
      double* values_diagonal = values + (worhp_w_.HM.nnz-nx_);

      // Initialize diagonal to zero
      for (int r=0; r<nx_; ++r) {
        values_diagonal[r] = 0.;
      }

      // Upper triangular part of the Hessian (note CCS -> CRS format change)
      for (int c=0; c<nx_; ++c) {
        for (int el=colind[c]; el<colind[c+1]; ++el) {
          if (row[el]>c) {
            // Strictly upper triangular
            *values_upper++ = data[el];
          } else if (row[el]==c) {
            // Diagonal separately
            values_diagonal[c] = data[el];
          }
        }
      }

      if (monitored("eval_h")) {
        std::cout << "x = " <<  hessLag.input(HESSLAG_X) << std::endl;
        std::cout << "obj_factor= " << obj_factor << std::endl;
        std::cout << "lambda = " << hessLag.input(HESSLAG_LAM_G) << std::endl;
        std::cout << "H = " << hessLag.output(HESSLAG_HESS) << std::endl;
      }

      if (regularity_check_ && !isRegular(hessLag.output(HESSLAG_HESS).data()))
          casadi_error("WorhpInterface::eval_h: NaN or Inf detected.");

      double time2 = clock();
      t_eval_h_ += (time2-time1)/CLOCKS_PER_SEC;
      n_eval_h_ += 1;
      log("eval_h ok");
      return true;
    } catch(exception& ex) {
      cerr << "eval_h failed: " << ex.what() << endl;
      return false;
    }
  }

  bool WorhpInterface::eval_jac_g(const double* x, double* values) {
    try {
      log("eval_jac_g started");

      // Quich finish if no constraints
      if (worhp_o_.m==0) {
        log("eval_jac_g quick return (m==0)");
        return true;
      }

      // Make sure generated
      casadi_assert(!jacG_.isNull());

      // Get Jacobian
      Function& jacG = this->jacG();

      double time1 = clock();

      // Pass the argument to the function
      jacG.setInput(x, JACG_X);
      jacG.setInput(input(NLP_SOLVER_P), JACG_P);

      // Evaluate the function
      jacG.evaluate();

      const DMatrix& J = jacG.output(JACG_JAC);

      std::copy(J.data().begin(), J.data().end(), values);

      if (monitored("eval_jac_g")) {
        cout << "x = " << jacG_.input().data() << endl;
        cout << "J = " << endl;
        jacG_.output().printSparse();
      }

      double time2 = clock();
      t_eval_jac_g_ += (time2-time1)/CLOCKS_PER_SEC;
      n_eval_jac_g_ += 1;
      log("eval_jac_g ok");
      return true;
    } catch(exception& ex) {
      cerr << "eval_jac_g failed: " << ex.what() << endl;
      return false;
    }
  }

  bool WorhpInterface::eval_f(const double* x, double scale, double& obj_value) {
    try {
      log("eval_f started");

      // Log time
      double time1 = clock();

      // Pass the argument to the function
      nlp_.setInput(x, NL_X);
      nlp_.setInput(input(NLP_SOLVER_P), NL_P);

      // Evaluate the function
      nlp_.evaluate();

      // Get the result
      nlp_.getOutput(obj_value, NL_F);

      // Printing
      if (monitored("eval_f")) {
        cout << "x = " << nlp_.input(NL_X) << endl;
        cout << "obj_value = " << obj_value << endl;
      }
      obj_value *= scale;

      if (regularity_check_ && !isRegular(nlp_.output().data()))
          casadi_error("WorhpInterface::eval_f: NaN or Inf detected.");

      double time2 = clock();
      t_eval_f_ += (time2-time1)/CLOCKS_PER_SEC;
      n_eval_f_ += 1;
      log("eval_f ok");
      return true;
    } catch(exception& ex) {
      cerr << "eval_f failed: " << ex.what() << endl;
      return false;
    }
  }

  bool WorhpInterface::eval_g(const double* x, double* g) {
    try {
      log("eval_g started");
      double time1 = clock();

      if (worhp_o_.m>0) {
        // Pass the argument to the function
        nlp_.setInput(x, NL_X);
        nlp_.setInput(input(NLP_SOLVER_P), NL_P);

        // Evaluate the function and tape
        nlp_.evaluate();

        // Ge the result
        nlp_.getOutput(g, NL_G);

        // Printing
        if (monitored("eval_g")) {
          cout << "x = " << nlp_.input(NL_X) << endl;
          cout << "g = " << nlp_.output(NL_G) << endl;
        }
      }

      if (regularity_check_ && !isRegular(nlp_.output(NL_G).data()))
          casadi_error("WorhpInterface::eval_g: NaN or Inf detected.");

      double time2 = clock();
      t_eval_g_ += (time2-time1)/CLOCKS_PER_SEC;
      n_eval_g_ += 1;
      log("eval_g ok");
      return true;
    } catch(exception& ex) {
      cerr << "eval_g failed: " << ex.what() << endl;
      return false;
    }
  }

  bool WorhpInterface::eval_grad_f(const double* x, double scale , double* grad_f) {
    try {
      log("eval_grad_f started");
      double time1 = clock();

      // Pass the argument to the function
      gradF_.setInput(x, NL_X);
      gradF_.setInput(input(NLP_SOLVER_P), NL_P);

      // Evaluate, adjoint mode
      gradF_.evaluate();

      // Get the result
      gradF_.output().get(grad_f, DENSE);

      // Scale
      for (int i=0; i<nx_; ++i) {
        grad_f[i] *= scale;
      }

      // Printing
      if (monitored("eval_grad_f")) {
        cout << "grad_f = " << gradF_.output() << endl;
      }

      if (regularity_check_ && !isRegular(gradF_.output().data()))
          casadi_error("WorhpInterface::eval_grad_f: NaN or Inf detected.");

      double time2 = clock();
      t_eval_grad_f_ += (time2-time1)/CLOCKS_PER_SEC;
      n_eval_grad_f_ += 1;
      // Check the result for regularity
      for (int i=0; i<nx_; ++i) {
        if (isnan(grad_f[i]) || isinf(grad_f[i])) {
          log("eval_grad_f: result not regular");
          return false;
        }
      }

      log("eval_grad_f ok");
      return true;
    } catch(exception& ex) {
      cerr << "eval_jac_f failed: " << ex.what() << endl;
      return false;
    }
  }

  void WorhpInterface::setOptionsFromFile(const std::string & file) {
    int status;
    char *cpy = new char[file.size()+1] ;
    strcpy(cpy, file.c_str());
    worhp_p_.initialised = true;
    ReadParamsNoInit(&status, cpy, &worhp_p_);
    delete cpy;

    for (int i=0;i<WorhpGetParamCount();++i) {
      WorhpType type = WorhpGetParamType(i+1);
      const char* name = WorhpGetParamName(i+1);
      if (strcmp(name, "Ares")==0) continue;
      switch (type) {
        case WORHP_BOOL_T:
          bool default_bool;
          WorhpGetBoolParam(&worhp_p_, name, &default_bool);
          setOption(WorhpGetParamName(i+1), default_bool);
          break;
        case WORHP_DOUBLE_T:
          double default_double;
          WorhpGetDoubleParam(&worhp_p_, name, &default_double);
          setOption(WorhpGetParamName(i+1), default_double);
          break;
        case WORHP_INT_T:
          int default_int;
          WorhpGetIntParam(&worhp_p_, name, &default_int);
          setOption(WorhpGetParamName(i+1), default_int);
          break;
        default:
          break; // do nothing
      }
    }


    setOption("qp_ipBarrier", worhp_p_.qp.ipBarrier);
    setOption("qp_ipComTol", worhp_p_.qp.ipComTol);
    setOption("qp_ipFracBound", worhp_p_.qp.ipFracBound);
    setOptionByEnumValue("qp_ipLsMethod", worhp_p_.qp.ipLsMethod);
    setOption("qp_ipMinAlpha", worhp_p_.qp.ipMinAlpha);
    setOption("qp_ipTryRelax", worhp_p_.qp.ipTryRelax);
    setOption("qp_ipRelaxDiv", worhp_p_.qp.ipRelaxDiv);
    setOption("qp_ipRelaxMult", worhp_p_.qp.ipRelaxMult);
    setOption("qp_ipRelaxMax", worhp_p_.qp.ipRelaxMax);
    setOption("qp_ipRelaxMin", worhp_p_.qp.ipRelaxMin);
    setOption("qp_ipResTol", worhp_p_.qp.ipResTol);
    setOption("qp_lsItMaxIter", worhp_p_.qp.lsItMaxIter);
    setOptionByEnumValue("qp_lsItMethod", worhp_p_.qp.lsItMethod);
    setOptionByEnumValue("qp_lsItPrecondMethod", worhp_p_.qp.lsItPrecondMethod);
    setOption("qp_lsRefineMaxIter", worhp_p_.qp.lsRefineMaxIter);
    setOption("qp_lsScale", worhp_p_.qp.lsScale);
    setOption("qp_lsTrySimple", worhp_p_.qp.lsTrySimple);
    setOption("qp_lsTol", worhp_p_.qp.lsTol);
    setOption("qp_maxIter", worhp_p_.qp.maxIter);
    setOptionByEnumValue("qp_method", worhp_p_.qp.method);
    setOption("qp_nsnBeta", worhp_p_.qp.nsnBeta);
    setOption("qp_nsnGradStep", worhp_p_.qp.nsnGradStep);
    setOption("qp_nsnKKT", worhp_p_.qp.nsnKKT);
    setOptionByEnumValue("qp_nsnLsMethod", worhp_p_.qp.nsnLsMethod);
    setOption("qp_nsnMinAlpha", worhp_p_.qp.nsnMinAlpha);
    setOption("qp_nsnSigma", worhp_p_.qp.nsnSigma);
    setOptionByEnumValue("qp_printLevel", worhp_p_.qp.printLevel);
    setOption("qp_scaleIntern", worhp_p_.qp.scaleIntern);
    setOption("qp_strict", worhp_p_.qp.strict);

    std::cout << "readparams status: " << status << std::endl;
  }

  map<int, string> WorhpInterface::calc_flagmap() {
  map<int, string> f;
  f[TerminateSuccess] = "TerminateSuccess";
  f[OptimalSolution] = "OptimalSolution";
  f[SearchDirectionZero] = "SearchDirectionZero";
  f[SearchDirectionSmall] = "SearchDirectionSmall";
  f[StationaryPointFound] = "StationaryPointFound";
  f[StationaryPointFound] = "StationaryPointFound";
  f[AcceptablePrevious] = "AcceptablePrevious";
  f[FritzJohn] = "FritzJohn";
  f[NotDiffable] = "NotDiffable";
  f[Unbounded] = "Unbounded";
  f[FeasibleSolution] = "FeasibleSolution";
  f[LowPassFilterOptimal] = "LowPassFilterOptimal";
  f[LowPassFilterAcceptable] = "LowPassFilterAcceptable";

  f[TerminateError] = "TerminateError";
  f[InitError] = "InitError";
  f[DataError] = "DataError";
  f[MaxCalls] = "MaxCalls";
  f[MaxIter] = "MaxIter";
  f[MinimumStepsize] = "MinimumStepsize";
  f[QPerror] = "QPerror";
  f[ProblemInfeasible] = "ProblemInfeasible";
  f[GroupsComposition] = "GroupsComposition";
  f[TooBig] = "TooBig";
  f[Timeout] = "Timeout";
  f[FDError] = "FDError";
  f[LocalInfeas] = "LocalInfeas";
  f[LicenseError] = "LicenseError";
  f[TerminatedByUser] = "TerminatedByUser";
  f[FunctionErrorF] = "FunctionErrorF";
  f[FunctionErrorG] = "FunctionErrorG";
  f[FunctionErrorDF] = "FunctionErrorDF";
  f[FunctionErrorDG] = "FunctionErrorDG";
  f[FunctionErrorHM] = "FunctionErrorHM";
  return f;
}

map<int, string> WorhpInterface::flagmap = WorhpInterface::calc_flagmap();

} // namespace casadi
