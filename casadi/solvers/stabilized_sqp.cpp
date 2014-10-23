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


#include "stabilized_sqp.hpp"
#include "casadi/core/std_vector_tools.hpp"
#include "casadi/core/matrix/sparsity_tools.hpp"
#include "casadi/core/matrix/matrix_tools.hpp"
#include "casadi/core/function/sx_function.hpp"
#include "casadi/core/sx/sx_tools.hpp"
#include "casadi/core/casadi_calculus.hpp"
#include <ctime>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cfloat>

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_NLPSOLVER_STABILIZEDSQP_EXPORT
  casadi_register_nlpsolver_stabilizedsqp(NlpSolverInternal::Plugin* plugin) {
    plugin->creator = StabilizedSqp::creator;
    plugin->name = "stabilizedsqp";
    plugin->doc = StabilizedSqp::meta_doc.c_str();
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOLVER_STABILIZEDSQP_EXPORT casadi_load_nlpsolver_stabilizedsqp() {
    NlpSolverInternal::registerPlugin(casadi_register_nlpsolver_stabilizedsqp);
  }

  StabilizedSqp::StabilizedSqp(const Function& nlp) : NlpSolverInternal(nlp) {
    casadi_warning("The SQP method is under development");
    addOption("stabilized_qp_solver",         OT_STRING,   GenericType(),
              "The Stabilized QP solver to be used by the SQP method");
    addOption("stabilized_qp_solver_options", OT_DICTIONARY, GenericType(),
              "Options to be passed to the Stabilized QP solver");
    addOption("hessian_approximation", OT_STRING, "exact",
              "limited-memory|exact");
    addOption("max_iter",           OT_INTEGER,     100,
              "Maximum number of SQP iterations");
    addOption("max_iter_ls",        OT_INTEGER,      20,
              "Maximum number of linesearch iterations");
    addOption("max_time",          OT_REAL,       1e12,
              "Timeout");
    addOption("tol_pr",            OT_REAL,       1e-5,
              "Stopping criterion for primal infeasibility");
    addOption("tol_du",            OT_REAL,       1e-5,
              "Stopping criterion for dual infeasability");
    addOption("c1",                OT_REAL,       0.001,
              "Armijo condition, coefficient of decrease in merit");
    addOption("beta",              OT_REAL,       0.5,
              "Line-search parameter, restoration factor of stepsize");
    addOption("merit_memory",      OT_INTEGER,      4,
              "Size of memory to store history of merit function values");
    addOption("lbfgs_memory",      OT_INTEGER,     10,
              "Size of L-BFGS memory.");
    addOption("regularize",        OT_BOOLEAN,  false,
              "Automatic regularization of Lagrange Hessian.");
    addOption("print_header",      OT_BOOLEAN,   true,
              "Print the header with problem statistics");
    addOption("min_step_size",     OT_REAL,   1e-10,
              "The size (inf-norm) of the step size should not become smaller than this.");

    // Monitors
    addOption("monitor",      OT_STRINGVECTOR, GenericType(),  "",
              "eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h|qp|dx", true);

    // New
    addOption("eps_active",        OT_REAL,      1e-6,
              "Threshold for the epsilon-active set.");
    addOption("nu",                OT_REAL,      1,
              "Parameter for primal-dual augmented Lagrangian.");
    addOption("phiWeight",         OT_REAL,      1e-5,
              "Weight used in pseudo-filter.");
    addOption("dvMax0",            OT_REAL,      100,
              "Parameter used to defined the max step length.");
    addOption("tau0",              OT_REAL,      1e-2,
              "Initial parameter for the merit function optimality threshold.");
    addOption("yEinitial",         OT_STRING,    "simple",
              "Initial multiplier. Simple (all zero) or least (LSQ).");
    addOption("alphaMin",          OT_REAL,      1e-3,
              "Used to check whether to increase rho.");
    addOption("sigmaMax",            OT_REAL,    1e+14,
              "Maximum rho allowed.");
    addOption("muR0",              OT_REAL,      1e-4,
              "Initial choice of regularization parameter");


    addOption("TReta1",            OT_REAL,      0.8,
              "Required predicted / actual decrease for TR increase");
    addOption("TReta2",            OT_REAL,      0.2,
              "Required predicted / actual decrease for TR decrease");
    addOption("gamma1",            OT_REAL,      2.,
              "Trust region increase parameter");
    addOption("gamma2",            OT_REAL,      1.,
              "Trust region update parameter");
    addOption("gamma3",            OT_REAL,      1.,
              "Trust region decrease parameter");

  }


  StabilizedSqp::~StabilizedSqp() {
  }

  void StabilizedSqp::init() {
    // Call the init method of the base class
    NlpSolverInternal::init();

    // Read options
    max_iter_ = getOption("max_iter");
    max_iter_ls_ = getOption("max_iter_ls");
    c1_ = getOption("c1");
    beta_ = getOption("beta");
    merit_memsize_ = getOption("merit_memory");
    lbfgs_memory_ = getOption("lbfgs_memory");
    tol_pr_ = getOption("tol_pr");
    tol_du_ = getOption("tol_du");
    regularize_ = getOption("regularize");
    exact_hessian_ = getOption("hessian_approximation")=="exact";
    min_step_size_ = getOption("min_step_size");

    eps_active_ = getOption("eps_active");
    nu_ = getOption("nu");
    phiWeight_ = getOption("phiWeight");
    dvMax_ = getOption("dvMax0");
    tau_ = getOption("tau0");
    alphaMin_ = getOption("alphaMin");
    sigmaMax_ = getOption("sigmaMax");
    muR_ = getOption("muR0");

    TReta1_ = getOption("TReta1");
    TReta2_ = getOption("TReta2");
    gamma1_ = getOption("gamma1");
    gamma2_ = getOption("gamma2");
    gamma3_ = getOption("gamma3");

    // Get/generate required functions
    gradF();
    jacG();
    if (exact_hessian_) {
      hessLag();
    }

    // Allocate a QP solver
    Sparsity H_sparsity = exact_hessian_ ? hessLag().output().sparsity()
        : Sparsity::dense(nx_, nx_);
    H_sparsity = H_sparsity + Sparsity::diag(nx_);
    Sparsity A_sparsity = jacG().isNull() ? Sparsity::sparse(0, nx_)
        : jacG().output().sparsity();

    std::string stabilized_qp_solver_name = getOption("stabilized_qp_solver");
    stabilized_qp_solver_ = StabilizedQpSolver(stabilized_qp_solver_name,
                                               qpStruct("h", H_sparsity, "a", A_sparsity));

    // Set options if provided
    if (hasSetOption("stabilized_qp_solver_options")) {
      Dictionary stabilized_qp_solver_options = getOption("stabilized_qp_solver_options");
      stabilized_qp_solver_.setOption(stabilized_qp_solver_options);
    }
    stabilized_qp_solver_.init();

    // Lagrange multipliers of the NLP
    mu_.resize(ng_);
    mu_cand_.resize(ng_);
    mu_x_.resize(nx_);
    mu_e_.resize(ng_);
    pi_.resize(ng_);
    pi2_.resize(ng_);

    // Lagrange gradient in the next iterate
    gLag_.resize(nx_);
    gLag_old_.resize(nx_);

    // Current linearization point
    x_.resize(nx_);
    x_cand_.resize(nx_);
    x_old_.resize(nx_);
    xtmp_.resize(nx_);

    // Constraint function value
    gk_.resize(ng_);
    gk_cand_.resize(ng_);
    gsk_.resize(ng_);
    gsk_cand_.resize(ng_);
    s_.resize(ng_);
    s_cand_.resize(ng_);

    // Hessian approximation
    Bk_ = DMatrix(H_sparsity);

    // Jacobian
    Jk_ = DMatrix(A_sparsity);


    // Bounds of the QP
    qp_LBA_.resize(ng_);
    qp_UBA_.resize(ng_);
    qp_LBX_.resize(nx_);
    qp_UBX_.resize(nx_);

    // QP solution
    dx_.resize(nx_);
    qp_DUAL_X_.resize(nx_);
    qp_DUAL_A_.resize(ng_);

    ds_.resize(ng_);
    dy_.resize(ng_);
    dv_.resize(ng_+nx_);

    // Merit function vectors
    dualpen_.resize(ng_);
    gradm_.resize(ng_+nx_);
    gradms_.resize(ng_);

    // Gradient of the objective
    gf_.resize(nx_);
    QPgf_.resize(nx_+ng_);

    // Primal-dual variables
    v_.resize(nx_+ng_);

    // Create Hessian update function
    if (!exact_hessian_) {
      // Create expressions corresponding to Bk, x, x_old, gLag and gLag_old
      SX Bk = SX::sym("Bk", H_sparsity);
      SX x = SX::sym("x", input(NLP_SOLVER_X0).sparsity());
      SX x_old = SX::sym("x", x.sparsity());
      SX gLag = SX::sym("gLag", x.sparsity());
      SX gLag_old = SX::sym("gLag_old", x.sparsity());

      SX sk = x - x_old;
      SX yk = gLag - gLag_old;
      SX qk = mul(Bk, sk);

      // Calculating theta
      SX skBksk = inner_prod(sk, qk);
      SX omega = if_else(inner_prod(yk, sk) < 0.2 * inner_prod(sk, qk),
                               0.8 * skBksk / (skBksk - inner_prod(sk, yk)),
                               1);
      yk = omega * yk + (1 - omega) * qk;
      SX theta = 1. / inner_prod(sk, yk);
      SX phi = 1. / inner_prod(qk, sk);
      SX Bk_new = Bk + theta * mul(yk, yk.T()) - phi * mul(qk, qk.T());

      // Inputs of the BFGS update function
      vector<SX> bfgs_in(BFGS_NUM_IN);
      bfgs_in[BFGS_BK] = Bk;
      bfgs_in[BFGS_X] = x;
      bfgs_in[BFGS_X_OLD] = x_old;
      bfgs_in[BFGS_GLAG] = gLag;
      bfgs_in[BFGS_GLAG_OLD] = gLag_old;
      bfgs_ = SXFunction(bfgs_in, Bk_new);
      bfgs_.init();

      // Initial Hessian approximation
      B_init_ = DMatrix::eye(nx_);
    }

    // Header
    if (static_cast<bool>(getOption("print_header"))) {
      cout << "-------------------------------------------" << endl;
      cout << "This is casadi::StabilizedSQPMethod." << endl;
      if (exact_hessian_) {
        cout << "Using exact Hessian" << endl;
      } else {
        cout << "Using limited memory BFGS Hessian approximation" << endl;
      }
      cout << endl;
      cout << "Number of variables:                       " << setw(9) << nx_ << endl;
      cout << "Number of constraints:                     " << setw(9) << ng_ << endl;
      cout << "Number of nonzeros in constraint Jacobian: " << setw(9) << A_sparsity.size() << endl;
      cout << "Number of nonzeros in Lagrangian Hessian:  " << setw(9) << H_sparsity.size() << endl;
      cout << endl;
    }
  }

  void StabilizedSqp::evaluate() {
    if (inputs_check_) checkInputs();
    checkInitialBounds();

    // Get problem data
    const vector<double>& x_init = input(NLP_SOLVER_X0).data();
    const vector<double>& lbx = input(NLP_SOLVER_LBX).data();
    const vector<double>& ubx = input(NLP_SOLVER_UBX).data();
    const vector<double>& lbg = input(NLP_SOLVER_LBG).data();
    const vector<double>& ubg = input(NLP_SOLVER_UBG).data();

    // Set linearization point to initial guess
    copy(x_init.begin(), x_init.end(), x_.begin());
    // make initial point feasible w.r.t. bounds
    for (int k=0; k<x_.size(); k++) {
        if (lbx[k] > x_[k])
            x_[k] = lbx[k];
        else if (x_[k] > ubx[k])
            x_[k] = ubx[k];
    }

    // Initialize Lagrange multipliers of the NLP
    copy(input(NLP_SOLVER_LAM_G0).begin(), input(NLP_SOLVER_LAM_G0).end(), mu_.begin());
    copy(input(NLP_SOLVER_LAM_X0).begin(), input(NLP_SOLVER_LAM_X0).end(), mu_x_.begin());

    // Initial constraint Jacobian
    eval_jac_g(x_, gk_, Jk_);

    // Initial objective gradient
    eval_grad_f(x_, fk_, gf_);

    normgf_ = norm_2(gf_);

    // Initialize or reset the Hessian or Hessian approximation
    reg_ = 0;
    if (exact_hessian_) {
      // eval_h(x_, mu_, 1.0, Bk_);
      // For first iteration, do not use
    } else {
      reset_h();
    }

    // Evaluate the initial gradient of the Lagrangian
    copy(gf_.begin(), gf_.end(), gLag_.begin());
    if (ng_>0) DMatrix::mul_no_alloc_tn(Jk_, mu_, gLag_);
    // gLag += mu_x_;
    transform(gLag_.begin(), gLag_.end(), mu_x_.begin(), gLag_.begin(), plus<double>());

    // Number of SQP iterations
    int iter = 0;

    // Number of line-search iterations
    int ls_iter = 0;

    // Last linesearch successfull
    bool ls_success = true;

    // Reset
    merit_mem_.clear();
    sigma_ = 1.;

    for (int i=0;i<ng_;++i)
      s_[i] = std::max(lbg[i], std::min(gk_[i], ubg[i]));//std::min(gk_[i]+mu_e_[i]*muR_, ubg[i]));
    for (int i=0;i<ng_;++i)
      gsk_[i] = gk_[i]-s_[i];

    TRsuccess_ = 0;
    TRDelta_ = std::max(norm_inf(s_), std::max(norm_inf(x_), std::max(norm_inf(gf_), 1.)));
//std::max(std::max(norm_inf(gf_), std::max(norm_inf(x_), std::max(norm_inf(s_), 1.))));
    rhoap_ = 0;
    rhoap_mu_ = 0;

    // Default stepsize
    double t = 0;

    // MAIN OPTIMIZATION LOOP
    double initial_time = clock();
    while (true) {

      // Primal infeasability
      double pr_inf = primalInfeasibility(x_, lbx, ubx, gk_, lbg, ubg);

      // inf-norm of lagrange gradient
      double gLag_norminf = norm_inf(gLag_);

      // inf-norm of step
      double dx_norminf = norm_inf(dx_);

      char info = ' ';

      // 1-norm of lagrange gradient
      //double gLag_norm1 = norm_1(gLag_);

      // 1-norm of step
      double dx_norm1 = norm_1(dx_);

      // Print header occasionally
      if (iter % 10 == 0) printIteration(cout);

      // This log entry is here to avoid #822
      log("Checking Stopping criteria");

      // Call callback function if present
      if (!callback_.isNull()) {
        if (!output(NLP_SOLVER_F).isEmpty()) output(NLP_SOLVER_F).set(fk_);
        if (!output(NLP_SOLVER_X).isEmpty()) output(NLP_SOLVER_X).set(x_);
        if (!output(NLP_SOLVER_LAM_G).isEmpty()) output(NLP_SOLVER_LAM_G).set(mu_);
        if (!output(NLP_SOLVER_LAM_X).isEmpty()) output(NLP_SOLVER_LAM_X).set(mu_x_);
        if (!output(NLP_SOLVER_G).isEmpty()) output(NLP_SOLVER_G).set(gk_);
        int ret = callback_(ref_, user_data_);

        if (!ret) {
          cout << endl;
          cout << "casadi::StabilizedSQPMethod: aborted by callback..." << endl;
          stats_["return_status"] = "User_Requested_Stop";
          break;
        }
      }
      normJ_ = norm1matrix(Jk_);
      // Default stepsize
      //double t = 0;



      normc_ = norm_2(gk_);
      normcs_ = norm_2(gsk_);

      scaleg_ = 1+normc_*normJ_;
      scaleglag_ = std::max(1., std::max(normgf_, (std::max(1., norm_2(mu_)) * normJ_)));

      if (exact_hessian_ && iter==0) {
        Bk_.setAll(0);
        Bk_(Sparsity::diag(nx_)) = 0.01 *scaleglag_;
      }

      // Checking convergence criteria
      if (pr_inf/scaleglag_ < tol_pr_ && gLag_norminf/scaleglag_ < tol_du_) {
        printIteration(cout, iter, fk_, pr_inf, gLag_norminf, dx_norm1, reg_,
                       TRDelta_, ls_iter, ls_success, ' ');
        cout << endl;
        cout << "casadi::StabilizedSQPMethod: Convergence achieved after "
             << iter << " iterations." << endl;
        stats_["return_status"] = "Solve_Succeeded";
        break;
      }

      if (iter==0) {
        phiMaxO_ = std::max(gLag_norminf+pr_inf+10., 1000.);
        phiMaxV_ = phiMaxO_;
        transform(mu_.begin(), mu_.end(), mu_e_.begin(), dualpen_.begin(), minus<double>());
        for (int i=0;i<ng_;++i)
        dualpen_[i] = -dualpen_[i]*(1/sigma_);
        transform(gsk_.begin(), gsk_.end(), dualpen_.begin(), dualpen_.begin(), plus<double>());
        Merit_ =
          fk_+inner_prod(mu_e_, gsk_)+.5*sigma_*(normcs_*normcs_+nu_*std::pow(norm_2(dualpen_), 2));

      } else { //do updates
        pr_inf = pr_inf / scaleglag_;
        gLag_norminf = gLag_norminf/scaleglag_;
        phiV_ = gLag_norminf+phiWeight_*pr_inf;
        phiO_ = phiWeight_*gLag_norminf+pr_inf;

        if (phiV_ <= 0.9*phiMaxV_) {
          phiMaxV_ = std::min(0.9*phiMaxV_, 0.75*phiV_+.25*phiMaxV_);
          info = 'O';
          copy(mu_.begin(), mu_.end(), mu_e_.begin());
          TRDelta_ = std::max(norm_inf(dx_)*3., TRDelta_);
        } else if (phiO_<=0.9*phiMaxO_) {
          phiMaxO_ = std::min(0.9*phiMaxO_, 0.75*phiO_+.25*phiMaxO_);
          info = 'V';
          copy(mu_.begin(), mu_.end(), mu_e_.begin());
          TRDelta_ = std::max(norm_inf(dx_)*3., TRDelta_);
        } else {
          transform(mu_x_.begin(), mu_x_.end(), gradm_.begin(), gradm_.begin(), plus<double>());

          double Mopt = std::max(norm_inf(gradm_), norm_inf(gradms_));

          if (Mopt<=tau_) {
            for (int i=0;i<ng_;++i)
            mu_e_[i] = std::max(-ymax, std::min(ymax, mu_[i]));
            info = 'M';
            tau_ = tau_/2;
            muR_ = std::min(muR_/2, std::pow(std::max(pr_inf, gLag_norminf), 0.8));
          } else {
            info = 'F';
          }
        }
        muR_ = std::min(muR_, std::pow(std::max(pr_inf, gLag_norminf), 0.8));

      }

      if (iter >= max_iter_) {
        cout << endl;
        cout << "casadi::StabilizedSQPMethod: Maximum number of iterations reached." << endl;
        stats_["return_status"] = "Maximum_Iterations_Exceeded";
        break;
      }

      if (iter > 0 && dx_norminf <= min_step_size_) {
        cout << endl;
        cout << "casadi::StabilizedSQPMethod: Search direction becomes too small without "
            "convergence criteria being met." << endl;
        stats_["return_status"] = "Search_Direction_Becomes_Too_Small";
        break;
      }

      if ((clock()-initial_time)/CLOCKS_PER_SEC > static_cast<double>(getOption("max_time"))) {
        cout << endl;
        cout << "casadi::StabilizedSQPMethod: Maximum time (" << getOption("max_time")
             << " sec.) exceeded." << endl;
        stats_["return_status"] = "Maximum_Time_Exceeded";
        break;
      }

      printIteration(cout, iter, fk_, pr_inf, gLag_norminf, dx_norm1, reg_,
                     TRDelta_, ls_iter, ls_success, info);

      // Start a new iteration
      iter++;

      log("Formulating QP");
      // Formulate the QP
      transform(lbx.begin(), lbx.end(), x_.begin(), qp_LBX_.begin(), minus<double>());
      transform(ubx.begin(), ubx.end(), x_.begin(), qp_UBX_.begin(), minus<double>());
      transform(lbg.begin(), lbg.end(), gk_.begin(), qp_LBA_.begin(), minus<double>());
      transform(ubg.begin(), ubg.end(), gk_.begin(), qp_UBA_.begin(), minus<double>());

      for (int i=0;i<nx_;++i) {
        qp_LBX_[i] = std::max(qp_LBX_[i], -TRDelta_);
        qp_UBX_[i] = std::min(qp_UBX_[i], TRDelta_);
      }

      // Solve the QP
      solve_QP(Bk_, gf_, qp_LBX_, qp_UBX_, Jk_, qp_LBA_, qp_UBA_, dx_,
               qp_DUAL_X_, qp_DUAL_A_, muR_, mu_, mu_e_);
      log("QP solved");

      for (int i=0;i<ng_;++i) {
        dy_[i] = qp_DUAL_A_[i]-mu_[i];
      }

      // Detecting indefiniteness
      double gain = quad_form(dx_, Bk_);

      mat_vec(dx_, Jk_, ds_);
      for (int i=0;i<ng_;++i) {
        ds_[i] = ds_[i]+gsk_[i]-muR_*(qp_DUAL_A_[i]-mu_e_[i]);
      }

      double muhat;
      //make sure, if nu=0 (so using classical augLag) that muR is small enough
      if (nu_==0) {
        mat_vectran(mu_e_, Jk_, xtmp_);
        transform(xtmp_.begin(), xtmp_.end(), gk_.begin(), xtmp_.begin(), plus<double>());
        muhat = inner_prod(xtmp_, dx_)-inner_prod(mu_e_, ds_)+.5*gain;
        muhat = inner_prod(gsk_, gsk_)/abs(muhat);
        transform(qp_DUAL_A_.begin(), qp_DUAL_A_.end(), mu_e_.begin(), pi2_.begin(),
                  minus<double>());
        muhat = std::min(muhat, norm_2(gsk_)/(2*norm_2(pi2_)));
        if (muR_>muhat)
        muR_ = muhat;
      }

      // Calculate line search values
      meritfg();


      //now calculate grads and inner product with both
      double rhsmerit = 0;
      for (int i=0;i<nx_;++i)
        rhsmerit = rhsmerit+ dx_[i]*gradm_[i];

      for (int i=0;i<ng_;++i)
        rhsmerit = rhsmerit+(dy_[i])*gradm_[i+nx_]+ds_[i]*gradms_[i];

      if (nu_== 0 && rhsmerit > 0) {
        for (int i=0;i<ng_;++i)
          rhsmerit = rhsmerit+gsk_[i]*(dy_[i]);
      }

      // Reset line-search counter, success marker
      ls_iter = 0;
      ls_success = true;
             // Calculate candidate Merit function
      t = 1; // std::min(1.0, dvMax_/std::max(norm_inf(dx_), norm_inf(dy_)));

      double dvHMdv = gain;

      mat_vec(dx_, Jk_, pi_);
      dvHMdv += (1+nu_)/muR_*std::pow(norm_2(pi_), 2);
      dvHMdv -= 2*(1+nu_)/muR_*inner_prod(pi_, ds_);
      dvHMdv -= 2*nu_*inner_prod(pi_, dy_);
      dvHMdv += (1+nu_)/muR_*std::pow(norm_2(ds_), 2);
      dvHMdv+=2*nu_*inner_prod(ds_, dy_);
      dvHMdv -= nu_*muR_*std::pow(norm_2(dy_), 2);

      for (int i=0; i<nx_; ++i) {
        x_cand_[i] = x_[i] + t * dx_[i];
      }

      // Evaluating objective and constraints
      eval_f(x_cand_, fk_cand_);
      eval_g(x_cand_, gk_cand_);

      for (int i=0;i<ng_;++i) {
        mu_cand_[i] = mu_[i]+t*(dy_[i]);
        s_cand_[i] = s_[i]+t*ds_[i];//std::min(ubg[i], std::max(lbg[i], gk_cand_[i]));//
        gsk_cand_[i] = gk_cand_[i]-s_cand_[i];
      }

      // Calculating merit-function in candidate
      normc_cand_ = norm_2(gk_cand_);
      normcs_cand_ = norm_2(gsk_cand_);

      transform(mu_cand_.begin(), mu_cand_.end(), mu_e_.begin(), dualpen_.begin(), minus<double>());
      for (int i=0;i<ng_;++i)
        dualpen_[i] = -dualpen_[i]*(1/sigma_);
      transform(gsk_cand_.begin(), gsk_cand_.end(), dualpen_.begin(), dualpen_.begin(),
                plus<double>());
      Merit_cand_ = fk_cand_+inner_prod(mu_e_, gsk_cand_)+.5*sigma_*(normcs_cand_*normcs_cand_+nu_*std::pow(norm_2(dualpen_), 2));   // NOLINT(whitespace/line_length)

      transform(mu_.begin(), mu_.end(), mu_e_.begin(), dualpen_.begin(), minus<double>());
      for (int i=0;i<ng_;++i)
             dualpen_[i] = -dualpen_[i]*(muR_);
      transform(gsk_.begin(), gsk_.end(), dualpen_.begin(), dualpen_.begin(), plus<double>());
      Merit_mu_ = fk_+inner_prod(mu_e_, gsk_)+.5*(1/muR_)*(normcs_*normcs_+nu_*std::pow(norm_2(dualpen_), 2));   // NOLINT(whitespace/line_length)

      transform(mu_cand_.begin(), mu_cand_.end(), mu_e_.begin(), dualpen_.begin(), minus<double>());
      for (int i=0;i<ng_;++i)
             dualpen_[i] = -dualpen_[i]*(muR_);
      transform(gsk_cand_.begin(), gsk_cand_.end(), dualpen_.begin(), dualpen_.begin(),
                plus<double>());
      Merit_mu_cand_ = fk_cand_+inner_prod(mu_e_, gsk_cand_)+.5*(1/muR_)*(normcs_cand_*normcs_cand_+nu_*std::pow(norm_2(dualpen_), 2));   // NOLINT(whitespace/line_length)


      if (Merit_cand_ - Merit_ > 0) {
        rhoap_ = 0;
      } else {
        rhoap_ = (Merit_cand_-Merit_) / (rhsmerit+0.5*dvHMdv);
      }

      rhoap_mu_ =
        (Merit_mu_cand_ - Merit_mu_ > 0)? 0 : (Merit_mu_cand_ - Merit_mu_) / (rhsmerit+0.5*dvHMdv);

      if ((rhoap_ > TReta2_) || (rhoap_mu_ > TReta2_)) {
        if ((rhoap_ > TReta1_ || rhoap_mu_ > TReta1_) &&
           (std::max(norm_inf(dx_), norm_inf(ds_))>0.01*TRDelta_)) {
          TRsuccess_++;
          TRDelta_ = std::sqrt(std::pow(TRDelta_*std::pow(gamma1_, TRsuccess_), 2)+std::pow(std::max(norm_inf(dx_), norm_inf(ds_))*std::pow(gamma1_, TRsuccess_), 2));  // NOLINT(whitespace/line_length)
        } else if (rhoap_ < TReta1_ && rhoap_mu_ <TReta1_) {
          TRsuccess_ = 0;
          TRDelta_ = TRDelta_*gamma2_;
        }

      } else {

        // Line-search
        log("Starting line-search");

        casadi_assert_message(max_iter_ls_ > 0, "max line search iterations should be > 0");

        // Line-search loop
        while (true) {

          ls_iter++;
          // Calculating maximal merit function value so far
          if (Merit_cand_ <= Merit_ + c1_*t * rhsmerit) {

            // Accepting candidate
            log("Line-search completed, candidate accepted");
            break;
          }


          // Line-search not successful, but we accept it.
          //do mu merit comparison as per flexible penalty

          if (Merit_mu_cand_ <= Merit_mu_+c1_*t*rhsmerit) {
            sigma_ = std::min(std::min(2*sigma_, 1/muR_), sigmaMax_);
            break;
          }

          if (ls_iter == max_iter_ls_) {
            ls_success = false;
            log("Line-search completed, maximum number of iterations");
            break;
          }

          // Backtracking
          t = beta_ * t;

          for (int i=0; i<nx_; ++i) {
            x_cand_[i] = x_[i] + t * dx_[i];
          }

          // Evaluating objective and constraints
          eval_f(x_cand_, fk_cand_);
          eval_g(x_cand_, gk_cand_);

          for (int i=0;i<ng_;++i) {
            mu_cand_[i] = mu_[i]+t*(dy_[i]);
            s_cand_[i] = s_[i]+t*ds_[i];//std::min(ubg[i], std::max(lbg[i], gk_cand_[i]));//
            gsk_cand_[i] = gk_cand_[i]-s_cand_[i];
          }

          // Calculating merit-function in candidate

          normc_cand_ = norm_2(gk_cand_);
          normcs_cand_ = norm_2(gsk_cand_);

          transform(mu_cand_.begin(), mu_cand_.end(), mu_e_.begin(), dualpen_.begin(),
                    minus<double>());
          for (int i=0;i<ng_;++i)
                 dualpen_[i] = -dualpen_[i]*(1/sigma_);
          transform(gsk_cand_.begin(), gsk_cand_.end(), dualpen_.begin(),
                    dualpen_.begin(), plus<double>());
          Merit_cand_ = fk_cand_+inner_prod(mu_e_, gsk_cand_)+.5*sigma_*(normcs_cand_*normcs_cand_+nu_*std::pow(norm_2(dualpen_), 2));  // NOLINT(whitespace/line_length)

          transform(mu_cand_.begin(), mu_cand_.end(), mu_e_.begin(), dualpen_.begin(),
                    minus<double>());
          for (int i=0;i<ng_;++i)
                 dualpen_[i] = -dualpen_[i]*(muR_);
          transform(gsk_cand_.begin(), gsk_cand_.end(), dualpen_.begin(), dualpen_.begin(),
                    plus<double>());
          Merit_mu_cand_ = fk_cand_+inner_prod(mu_e_, gsk_cand_)+.5*(1/muR_)*(normcs_cand_*normcs_cand_+nu_*std::pow(norm_2(dualpen_), 2));  // NOLINT(whitespace/line_length)

        } // while (true)
        if (ls_success) {
          TRDelta_ = gamma3_*t*std::max(norm_inf(dx_), norm_inf(ds_));
        } else {
          TRDelta_ = TRDelta_/25.;
        }
      }

      // Candidate accepted, update dual variables
      for (int i=0; i<ng_; ++i) mu_[i] = t * dy_[i] + mu_[i];
      for (int i=0; i<nx_; ++i) mu_x_[i] = t * qp_DUAL_X_[i] + (1 - t) * mu_x_[i];
      for (int i=0;i<ng_;++i) s_[i] = t*ds_[i]+s_[i];

      if (!exact_hessian_) {
        // Evaluate the gradient of the Lagrangian with the old x but new mu (for BFGS)
        copy(gf_.begin(), gf_.end(), gLag_old_.begin());
        if (ng_>0) DMatrix::mul_no_alloc_tn(Jk_, mu_, gLag_old_);
        // gLag_old += mu_x_;
        transform(gLag_old_.begin(), gLag_old_.end(), mu_x_.begin(), gLag_old_.begin(),
                  plus<double>());
      }
      // Candidate accepted, update the primal variable
      copy(x_.begin(), x_.end(), x_old_.begin());
      copy(x_cand_.begin(), x_cand_.end(), x_.begin());
      Merit_ = Merit_cand_;

      // Evaluate the constraint Jacobian
      log("Evaluating jac_g");
      eval_jac_g(x_, gk_, Jk_);

      // Evaluate the gradient of the objective function
      log("Evaluating grad_f");
      eval_grad_f(x_, fk_, gf_);

      // Evaluate the gradient of the Lagrangian with the new x and new mu
      copy(gf_.begin(), gf_.end(), gLag_.begin());
      if (ng_>0) DMatrix::mul_no_alloc_tn(Jk_, mu_, gLag_);
      // gLag += mu_x_;
      transform(gLag_.begin(), gLag_.end(), mu_x_.begin(), gLag_.begin(), plus<double>());

      // Updating Lagrange Hessian
      if (!exact_hessian_) {
        log("Updating Hessian (BFGS)");
        // BFGS with careful updates and restarts
        if (iter % lbfgs_memory_ == 0) {
          // Reset Hessian approximation by dropping all off-diagonal entries
          const vector<int>& colind = Bk_.colind();      // Access sparsity (column offset)
          const vector<int>& row = Bk_.row();            // Access sparsity (row)
          vector<double>& data = Bk_.data();             // Access nonzero elements
          for (int cc=0; cc<colind.size()-1; ++cc) {       // Loop over the columns of the Hessian
            for (int el=colind[cc]; el<colind[cc+1]; ++el) {
              // Loop over the nonzero elements of the column
              if (cc!=row[el]) data[el] = 0;               // Remove if off-diagonal entries
            }
          }
        }

        // Pass to BFGS update function
        bfgs_.setInput(Bk_, BFGS_BK);
        bfgs_.setInput(x_, BFGS_X);
        bfgs_.setInput(x_old_, BFGS_X_OLD);
        bfgs_.setInput(gLag_, BFGS_GLAG);
        bfgs_.setInput(gLag_old_, BFGS_GLAG_OLD);

        // Update the Hessian approximation
        bfgs_.evaluate();

        // Get the updated Hessian
        bfgs_.getOutput(Bk_);
      } else {
        // Exact Hessian
        log("Evaluating hessian");
        eval_h(x_, mu_, 1.0, Bk_);
      }

      for (int i=0;i<ng_;++i)
        gsk_[i] = gk_[i]-s_[i];
      //normc_ = norm_2(gk_);
      //normcs_ = norm_2(gsk_);
      dvMax_ = std::max(std::min(std::pow(beta_, (ls_iter-1))*dvMax_, 100.), 1e-8);
    }

    // Save results to outputs
    output(NLP_SOLVER_F).set(fk_);
    output(NLP_SOLVER_X).set(x_);
    output(NLP_SOLVER_LAM_G).set(mu_);
    output(NLP_SOLVER_LAM_X).set(mu_x_);
    output(NLP_SOLVER_G).set(gk_);

    // Save statistics
    stats_["iter_count"] = iter;
  }

  void StabilizedSqp::printIteration(std::ostream &stream) {
    stream << setw(4)  << "iter";
    stream << setw(14) << "objective";
    stream << setw(9) << "inf_pr";
    stream << setw(9) << "inf_du";
    stream << setw(9) << "||d||";
    stream << setw(9) << "TRDelta";
    stream << setw(7) << "lg(rg)";
    stream << setw(3) << "ls";
    stream << setw(2) << "info";
    //stream << ' ';
    stream << endl;
  }

  void StabilizedSqp::printIteration(std::ostream &stream, int iter, double obj,
                                             double pr_inf, double du_inf,
                                             double dx_norm, double rg, double TRdelta,
                                             int ls_trials, bool ls_success, char info) {
    stream << setw(4) << iter;
    stream << scientific;
    stream << setw(14) << setprecision(6) << obj;
    stream << setw(9) << setprecision(2) << pr_inf;
    stream << setw(9) << setprecision(2) << du_inf;
    stream << setw(9) << setprecision(2) << dx_norm;
    stream << setw(9) << setprecision(2) << TRdelta;
    stream << fixed;
    if (rg>0) {
      stream << setw(7) << setprecision(2) << log10(rg);
    } else {
      stream << setw(7) << "-";
    }
    stream << setw(3) << ls_trials;
    stream << (ls_success ? ' ' : 'F');
    stream << setw(3)<< info;
    stream << endl;
  }

  double StabilizedSqp::quad_form(const std::vector<double>& x, const DMatrix& A) {
    // Assert dimensions
    casadi_assert(x.size()==A.size1() && x.size()==A.size2());

    // Access the internal data of A
    const std::vector<int> &A_colind = A.colind();
    const std::vector<int> &A_row = A.row();
    const std::vector<double> &A_data = A.data();

    // Return value
    double ret=0;

    // Loop over the columns of A
    for (int cc=0; cc<x.size(); ++cc) {
      // Loop over the nonzeros of A
      for (int el=A_colind[cc]; el<A_colind[cc+1]; ++el) {
        // Get column
        int rr = A_row[el];

        // Add contribution
        ret += x[cc]*A_data[el]*x[rr];
      }
    }

    return ret;
  }

  void StabilizedSqp::reset_h() {
    // Initial Hessian approximation of BFGS
    if (!exact_hessian_) {
      Bk_.set(B_init_);
    }

    if (monitored("eval_h")) {
      cout << "x = " << x_ << endl;
      cout << "H = " << endl;
      Bk_.printSparse();
    }
  }

  double StabilizedSqp::getRegularization(const Matrix<double>& H) {
    const vector<int>& colind = H.colind();
    const vector<int>& row = H.row();
    const vector<double>& data = H.data();
    double reg_param = 0;
    for (int cc=0; cc<colind.size()-1; ++cc) {
      double mineig = 0;
      for (int el=colind[cc]; el<colind[cc+1]; ++el) {
        int rr = row[el];
        if (rr == cc) {
          mineig += data[el];
        } else {
          mineig -= fabs(data[el]);
        }
      }
      reg_param = fmin(reg_param, mineig);
    }
    return -reg_param;
  }

  void StabilizedSqp::regularize(Matrix<double>& H, double reg) {
    const vector<int>& colind = H.colind();
    const vector<int>& row = H.row();
    vector<double>& data = H.data();

    for (int cc=0; cc<colind.size()-1; ++cc) {
      for (int el=colind[cc]; el<colind[cc+1]; ++el) {
        int rr = row[el];
        if (rr==cc) {
          data[el] += reg;
        }
      }
    }
  }


  void StabilizedSqp::eval_h(const std::vector<double>& x,
                                     const std::vector<double>& lambda, double sigma,
                                     Matrix<double>& H) {
    try {
      // Get function
      Function& hessLag = this->hessLag();

      // Pass the argument to the function
      hessLag.setInput(x, HESSLAG_X);
      hessLag.setInput(input(NLP_SOLVER_P), HESSLAG_P);
      hessLag.setInput(sigma, HESSLAG_LAM_F);
      hessLag.setInput(lambda, HESSLAG_LAM_G);

      // Evaluate
      hessLag.evaluate();

      // Get results
      hessLag.getOutput(H);

      if (monitored("eval_h")) {
        cout << "x = " << x << endl;
        cout << "H = " << endl;
        H.printSparse();
      }

      // Determing regularization parameter with Gershgorin theorem
      if (regularize_) {
        reg_ = getRegularization(H);
        if (reg_ > 0) {
          regularize(H, reg_);
        }
      }

    } catch(exception& ex) {
      cerr << "eval_h failed: " << ex.what() << endl;
      throw;
    }
  }

  void StabilizedSqp::eval_g(const std::vector<double>& x, std::vector<double>& g) {
    try {

      // Quick return if no constraints
      if (ng_==0) return;

      // Pass the argument to the function
      nlp_.setInput(x, NL_X);
      nlp_.setInput(input(NLP_SOLVER_P), NL_P);

      // Evaluate the function and tape
      nlp_.evaluate();

      // Ge the result
      nlp_.output(NL_G).get(g, DENSE);

      // Printing
      if (monitored("eval_g")) {
        cout << "x = " << nlp_.input(NL_X) << endl;
        cout << "g = " << nlp_.output(NL_G) << endl;
      }
    } catch(exception& ex) {
      cerr << "eval_g failed: " << ex.what() << endl;
      throw;
    }
  }

  void StabilizedSqp::eval_jac_g(const std::vector<double>& x,
                                         std::vector<double>& g, Matrix<double>& J) {
    try {
      // Quich finish if no constraints
      if (ng_==0) return;

      // Get function
      Function& jacG = this->jacG();

      // Pass the argument to the function
      jacG.setInput(x, NL_X);
      jacG.setInput(input(NLP_SOLVER_P), NL_P);

      // Evaluate the function
      jacG.evaluate();

      // Get the output
      jacG.output(1+NL_G).get(g, DENSE);
      jacG.output().get(J);

      if (monitored("eval_jac_g")) {
        cout << "x = " << x << endl;
        cout << "g = " << g << endl;
        cout << "J = " << endl;
        J.printSparse();
      }
    } catch(exception& ex) {
      cerr << "eval_jac_g failed: " << ex.what() << endl;
      throw;
    }
  }

  void StabilizedSqp::eval_grad_f(const std::vector<double>& x,
                                          double& f, std::vector<double>& grad_f) {
    try {
      // Get function
      Function& gradF = this->gradF();

      // Pass the argument to the function
      gradF.setInput(x, NL_X);
      gradF.setInput(input(NLP_SOLVER_P), NL_P);

      // Evaluate, adjoint mode
      gradF.evaluate();

      // Get the result
      gradF.output().get(grad_f, DENSE);
      gradF.output(1+NL_X).get(f);

      // Printing
      if (monitored("eval_f")) {
        cout << "x = " << x << endl;
        cout << "f = " << f << endl;
      }

      if (monitored("eval_grad_f")) {
        cout << "x      = " << x << endl;
        cout << "grad_f = " << grad_f << endl;
      }
    } catch(exception& ex) {
      cerr << "eval_grad_f failed: " << ex.what() << endl;
      throw;
    }
  }

  void StabilizedSqp::eval_f(const std::vector<double>& x, double& f) {
    try {
      // Pass the argument to the function
      nlp_.setInput(x, NL_X);
      nlp_.setInput(input(NLP_SOLVER_P), NL_P);

      // Evaluate the function
      nlp_.evaluate();

      // Get the result
      nlp_.getOutput(f, NL_F);

      // Printing
      if (monitored("eval_f")) {
        cout << "x = " << nlp_.input(NL_X) << endl;
        cout << "f = " << f << endl;
      }
    } catch(exception& ex) {
      cerr << "eval_f failed: " << ex.what() << endl;
      throw;
    }
  }

  void StabilizedSqp::solve_QP(const Matrix<double>& H, const std::vector<double>& g,
                                       const std::vector<double>& lbx,
                                       const std::vector<double>& ubx,
                                       const Matrix<double>& A,
                                       const std::vector<double>& lbA,
                                       const std::vector<double>& ubA,
                                       std::vector<double>& x_opt,
                                       std::vector<double>& lambda_x_opt,
                                       std::vector<double>& lambda_A_opt,
                                       double muR,
                                       const std::vector<double> & mu,
                                       const std::vector<double> & muE) {

    // Pass data to QP solver
    stabilized_qp_solver_.setInput(H, STABILIZED_QP_SOLVER_H);
    stabilized_qp_solver_.setInput(g, STABILIZED_QP_SOLVER_G);
    stabilized_qp_solver_.setInput(mu, STABILIZED_QP_SOLVER_MU);
    stabilized_qp_solver_.setInput(muE, STABILIZED_QP_SOLVER_MUE);
    stabilized_qp_solver_.setInput(muR, STABILIZED_QP_SOLVER_MUR);

    // Hot-starting if possible
    stabilized_qp_solver_.setInput(x_opt, STABILIZED_QP_SOLVER_X0);

    //TODO(Joel): Fix hot-starting of dual variables
    //qp_solver_.setInput(lambda_A_opt, QP_SOLVER_LAMBDA_INIT);

    // Pass simple bounds
    stabilized_qp_solver_.setInput(lbx, STABILIZED_QP_SOLVER_LBX);
    stabilized_qp_solver_.setInput(ubx, STABILIZED_QP_SOLVER_UBX);

    // Pass linear bounds
    if (ng_>0) {
      stabilized_qp_solver_.setInput(A, STABILIZED_QP_SOLVER_A);
      stabilized_qp_solver_.setInput(lbA, STABILIZED_QP_SOLVER_LBA);
      stabilized_qp_solver_.setInput(ubA, STABILIZED_QP_SOLVER_UBA);
    }

    if (monitored("qp")) {
      cout << "H = " << endl;
      H.printDense();
      cout << "A = " << endl;
      A.printDense();
      cout << "g = " << g << endl;
      cout << "lbx = " << lbx << endl;
      cout << "ubx = " << ubx << endl;
      cout << "lbA = " << lbA << endl;
      cout << "ubA = " << ubA << endl;
    }

    // Solve the QP
    stabilized_qp_solver_.evaluate();

    // Get the optimal solution
    stabilized_qp_solver_.getOutput(x_opt, QP_SOLVER_X);
    stabilized_qp_solver_.getOutput(lambda_x_opt, QP_SOLVER_LAM_X);
    stabilized_qp_solver_.getOutput(lambda_A_opt, QP_SOLVER_LAM_A);
    if (monitored("dx")) {
      cout << "dx = " << x_opt << endl;
    }

   }

  double StabilizedSqp::norm1matrix(const DMatrix& A) {
    // Access the arrays
    const std::vector<double>& v = A.data();
    const std::vector<int>& colind = A.colind();
    double ret = 0;
    std::vector<double> sums(A.size2(), 0);
    for (int cc=0; cc<colind.size()-1; ++cc) {
      double colsum = 0;
      for (int el=colind[cc]; el<colind[cc+1]; ++el) {
        colsum += abs(v[el]);
      }
      ret = max(ret, colsum);
    }
    return ret;
  }

  double StabilizedSqp::primalInfeasibility(const std::vector<double>& x,
                                                    const std::vector<double>& lbx,
                                                    const std::vector<double>& ubx,
                                                    const std::vector<double>& g,
                                                    const std::vector<double>& lbg,
                                                    const std::vector<double>& ubg) {
    // Linf-norm of the primal infeasibility
    double pr_inf = 0;

    // Bound constraints
    for (int j=0; j<x.size(); ++j) {
      pr_inf = max(pr_inf, lbx[j] - x[j]);
      pr_inf = max(pr_inf, x[j] - ubx[j]);
    }

    // Nonlinear constraints
    for (int j=0; j<g.size(); ++j) {
      pr_inf = max(pr_inf, lbg[j] - g[j]);
      pr_inf = max(pr_inf, g[j] - ubg[j]);
    }

    return pr_inf;
  }

void StabilizedSqp::mat_vectran(const std::vector<double>& x,
                                        const DMatrix& A, std::vector<double>& y) {
    // Assert dimensions
    casadi_assert(x.size()==A.size1() && y.size()==A.size2());

    // Access the internal data of A
    const std::vector<int> &A_colind = A.colind();
    const std::vector<int> &A_row = A.row();
    const std::vector<double> &A_data = A.data();


    for (int i=0;i<y.size();++i)
      y[i] = 0;
    // Loop over the columns of A
    for (int cc=0; cc<A.size2(); ++cc) {
      // Loop over the nonzeros of A
      for (int el=A_colind[cc]; el<A_colind[cc+1]; el++) {
        // Get row
        int rr = A_row[el];

        // Add contribution
        y[cc] += A_data[el]*x[rr];
      }
    }

  }

  void StabilizedSqp::mat_vec(const std::vector<double>& x,
                                      const DMatrix& A, std::vector<double>& y) {
    // Assert dimensions
    casadi_assert(x.size()==A.size2() && y.size()==A.size1());

    // Access the internal data of A
    const std::vector<int> &A_colind = A.colind();
    const std::vector<int> &A_row = A.row();
    const std::vector<double> &A_data = A.data();


    for (int i=0;i<y.size();++i)
      y[i] = 0;
    // Loop over the rows of A
    for (int cc=0; cc<A.size2(); ++cc) {
      // Loop over the nonzeros of A
      for (int el=A_colind[cc]; el<A_colind[cc+1]; el++) {
        // Get column
        int rr = A_row[el];

        // Add contribution
        y[rr] += A_data[el]*x[cc];
      }
    }

  }

  void StabilizedSqp::meritfg() {
    for (int i=0;i<ng_;++i) {
      pi_[i] = (1+nu_)/muR_*gsk_[i]+(1+nu_)*mu_e_[i]-nu_*mu_[i];
    }
    //transform(mu_e_.begin(), mu_e_.end(), pi_.begin(), pi_.begin(), plus<double>());
    //transform(pi_.begin(), pi_.end(), mu_.begin(), pi_.begin(), plus<double>());
    //for (int i=0;i<ng_;++i) {
    //  pi_[i] = nu_*pi_[i];
    //  pi2_[i] = gsk_[i]/muR_+(1-nu_)*mu_e_[i];
    //}
    //transform(pi_.begin(), pi_.end(), mu_e_.begin(), pi_.begin(), plus<double>());
    //transform(pi_.begin(), pi_.end(), pi2_.begin(), pi_.begin(), plus<double>());
    copy(gf_.begin(), gf_.end(), gradm_.begin());
    mat_vectran(pi_, Jk_, xtmp_);
    transform(xtmp_.begin(), xtmp_.end(), gradm_.begin(), gradm_.begin(), plus<double>());
    for (int i=0;i<ng_;++i) {
      gradm_[nx_+i] = nu_*(gsk_[i]-muR_*(mu_[i]-mu_e_[i]));
      gradms_[i] = -(mu_e_[i]+(1+nu_)*gsk_[i]/muR_-nu_*(mu_[i]-mu_e_[i]));
    }


  }

} // namespace casadi
