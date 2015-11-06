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


#include "sqpmethod.hpp"

#include "casadi/core/std_vector_tools.hpp"
#include "casadi/core/casadi_calculus.hpp"

#include <ctime>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cfloat>

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_NLPSOL_SQPMETHOD_EXPORT
      casadi_register_nlpsol_sqpmethod(Nlpsol::Plugin* plugin) {
    plugin->creator = Sqpmethod::creator;
    plugin->name = "sqpmethod";
    plugin->doc = Sqpmethod::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOL_SQPMETHOD_EXPORT casadi_load_nlpsol_sqpmethod() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_sqpmethod);
  }

  Sqpmethod::Sqpmethod(const std::string& name, const XProblem& nlp)
    : Nlpsol(name, nlp) {

    casadi_warning("The SQP method is under development");
    addOption("qpsol",         OT_STRING,   GenericType(),
              "The QP solver to be used by the SQP method");
    addOption("qpsol_options", OT_DICT, GenericType(),
              "Options to be passed to the QP solver");
    addOption("hessian_approximation", OT_STRING, "exact",
              "limited-memory|exact");
    addOption("max_iter",           OT_INTEGER,      50,
              "Maximum number of SQP iterations");
    addOption("max_iter_ls",        OT_INTEGER,       3,
              "Maximum number of linesearch iterations");
    addOption("tol_pr",            OT_REAL,       1e-6,
              "Stopping criterion for primal infeasibility");
    addOption("tol_du",            OT_REAL,       1e-6,
              "Stopping criterion for dual infeasability");
    addOption("c1",                OT_REAL,       1E-4,
              "Armijo condition, coefficient of decrease in merit");
    addOption("beta",              OT_REAL,       0.8,
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
              "eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h|qp|dx|bfgs", true);
    addOption("print_time",         OT_BOOLEAN,       true,
              "Print information about execution time");
  }


  Sqpmethod::~Sqpmethod() {
  }

  void Sqpmethod::init() {
    // Call the init method of the base class
    Nlpsol::init();

    // Read options
    max_iter_ = option("max_iter");
    max_iter_ls_ = option("max_iter_ls");
    c1_ = option("c1");
    beta_ = option("beta");
    merit_memsize_ = option("merit_memory");
    lbfgs_memory_ = option("lbfgs_memory");
    tol_pr_ = option("tol_pr");
    tol_du_ = option("tol_du");
    regularize_ = option("regularize");
    exact_hessian_ = option("hessian_approximation")=="exact";
    min_step_size_ = option("min_step_size");

    // Get/generate required functions
    gradF();
    jacG();
    if (exact_hessian_) {
      hessLag();
    }

    // Allocate a QP solver
    Hsp_ = exact_hessian_ ? hessLag().sparsity_out(0) : Sparsity::dense(nx_, nx_);
    Hsp_ = Hsp_ + Sparsity::diag(nx_);
    Asp_ = jacG().isNull() ? Sparsity(0, nx_) : jacG().sparsity_out(0);

    // QP solver options
    Dict qpsol_options;
    if (hasSetOption("qpsol_options")) {
      qpsol_options = option("qpsol_options");
    }

    // Allocate a QP solver
    qpsol_ = Function::qpsol("qpsol", option("qpsol"),
                             {{"h", Hsp_}, {"a", Asp_}},
                             qpsol_options);
    alloc(qpsol_);

    // Create Hessian update function
    if (!exact_hessian_) {
      // Create expressions corresponding to Bk, x, x_old, gLag and gLag_old
      SX Bk = SX::sym("Bk", Hsp_);
      SX x = SX::sym("x", sparsity_in(NLPSOL_X0));
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
      bfgs_ = Function("bfgs", bfgs_in, {Bk_new});
      alloc(bfgs_);

      // Initial Hessian approximation
      B_init_ = DMatrix(Hsp_);
      B_init_.set(DMatrix::eye(nx_));
    }

    // Header
    if (static_cast<bool>(option("print_header"))) {
      userOut()
        << "-------------------------------------------" << endl
        << "This is casadi::SQPMethod." << endl;
      if (exact_hessian_) {
        userOut() << "Using exact Hessian" << endl;
      } else {
        userOut() << "Using limited memory BFGS Hessian approximation" << endl;
      }
      userOut()
        << endl
        << "Number of variables:                       " << setw(9) << nx_ << endl
        << "Number of constraints:                     " << setw(9) << ng_ << endl
        << "Number of nonzeros in constraint Jacobian: " << setw(9) << Asp_.nnz() << endl
        << "Number of nonzeros in Lagrangian Hessian:  " << setw(9) << Hsp_.nnz() << endl
        << endl;
    }

    // Lagrange multipliers of the NLP
    alloc_w(ng_, true); // mu_
    alloc_w(nx_, true); // mu_x_

    // Current linearization point
    alloc_w(nx_, true); // x_
    alloc_w(nx_, true); // x_cand_
    alloc_w(nx_, true); // x_old_

    // Lagrange gradient in the next iterate
    alloc_w(nx_, true); // gLag_
    alloc_w(nx_, true); // gLag_old_

    // Constraint function value
    alloc_w(ng_, true); // gk_
    alloc_w(ng_, true); // gk_cand_

    // Gradient of the objective
    alloc_w(nx_, true); // gf_

    // Bounds of the QP
    alloc_w(ng_, true); // qp_LBA
    alloc_w(ng_, true); // qp_UBA
    alloc_w(nx_, true); // qp_LBX
    alloc_w(nx_, true); // qp_UBX

    // QP solution
    alloc_w(nx_, true); // dx_
    alloc_w(nx_, true); // qp_DUAL_X_
    alloc_w(ng_, true); // qp_DUAL_A_

    // Hessian approximation
    alloc_w(Hsp_.nnz(), true); // Bk_

    // Jacobian
    alloc_w(Asp_.nnz(), true); // Jk_
  }

  void Sqpmethod::evalD(void* mem, const double** arg, double** res, int* iw, double* w) {
    // Get input pointers
    x0_ = arg[NLPSOL_X0];
    p_ = arg[NLPSOL_P];
    lbx_ = arg[NLPSOL_LBX];
    ubx_ = arg[NLPSOL_UBX];
    lbg_ = arg[NLPSOL_LBG];
    ubg_ = arg[NLPSOL_UBG];
    lam_x0_ = arg[NLPSOL_LAM_X0];
    lam_g0_ = arg[NLPSOL_LAM_G0];

    // Get output pointers
//    x_ = res[NLPSOL_X];
    f_ = res[NLPSOL_F];
    g_ = res[NLPSOL_G];
    lam_x_ = res[NLPSOL_LAM_X];
    lam_g_ = res[NLPSOL_LAM_G];
    lam_p_ = res[NLPSOL_LAM_P];

    // Lagrange multipliers of the NLP
    mu_ = w; w += ng_;
    mu_x_ = w; w += nx_;

    // Current linearization point
    x_ = w; w += nx_;
    x_cand_ = w; w += nx_;
    x_old_ = w; w += nx_;

    // Lagrange gradient in the next iterate
    gLag_ = w; w += nx_;
    gLag_old_ = w; w += nx_;

    // Constraint function value
    gk_ = w; w += ng_;
    gk_cand_ = w; w += ng_;

    // Gradient of the objective
    gf_ = w; w += nx_;

    // Bounds of the QP
    qp_LBA_ = w; w += ng_;
    qp_UBA_ = w; w += ng_;
    qp_LBX_ = w; w += nx_;
    qp_UBX_ = w; w += nx_;

    // QP solution
    dx_ = w; w += nx_;
    qp_DUAL_X_ = w; w += nx_;
    qp_DUAL_A_ = w; w += ng_;

    // Hessian approximation
    Bk_ = w; w += Hsp_.nnz();

    // Jacobian
    Jk_ = w; w += Asp_.nnz();

    // Get work vectors
    arg1_ = arg + NLPSOL_NUM_IN;
    res1_ = res + NLPSOL_NUM_OUT;
    iw_ = iw;
    w_ = w;

    // Pass the inputs to the function
    for (int i=0; i<n_in(); ++i) {
      if (arg[i] != 0) {
        setInputNZ(arg[i], i);
      } else {
        setInput(0., i);
      }
    }

    if (inputs_check_) checkInputs();
    checkInitialBounds();

    if (gather_stats_) {
      Dict iterations;
      iterations["inf_pr"] = std::vector<double>();
      iterations["inf_du"] = std::vector<double>();
      iterations["ls_trials"] = std::vector<double>();
      iterations["d_norm"] = std::vector<double>();
      iterations["obj"] = std::vector<double>();
      stats_["iterations"] = iterations;
    }


    // Get problem data
    const vector<double>& x_init = input(NLPSOL_X0).data();
    const vector<double>& lbx = input(NLPSOL_LBX).data();
    const vector<double>& ubx = input(NLPSOL_UBX).data();
    const vector<double>& lbg = input(NLPSOL_LBG).data();
    const vector<double>& ubg = input(NLPSOL_UBG).data();

    // Set linearization point to initial guess
    copy(x_init.begin(), x_init.end(), x_);

    // Initialize Lagrange multipliers of the NLP
    copy(input(NLPSOL_LAM_G0)->begin(), input(NLPSOL_LAM_G0)->end(), mu_);
    copy(input(NLPSOL_LAM_X0)->begin(), input(NLPSOL_LAM_X0)->end(), mu_x_);

    t_eval_f_ = t_eval_grad_f_ = t_eval_g_ = t_eval_jac_g_ = t_eval_h_ =
        t_callback_fun_ = t_callback_prepare_ = t_mainloop_ = 0;

    n_eval_f_ = n_eval_grad_f_ = n_eval_g_ = n_eval_jac_g_ = n_eval_h_ = 0;

    double time1 = clock();

    // Initial constraint Jacobian
    eval_jac_g(x_, gk_, Jk_);

    // Initial objective gradient
    eval_grad_f(x_, &fk_, gf_);

    // Initialize or reset the Hessian or Hessian approximation
    reg_ = 0;
    if (exact_hessian_) {
      eval_h(x_, mu_, 1.0, Bk_);
    } else {
      reset_h();
    }

    // Evaluate the initial gradient of the Lagrangian
    copy(gf_, gf_+nx_, gLag_);
    if (ng_>0) casadi_mv_t(Jk_, Asp_, mu_, gLag_);
    // gLag += mu_x_;
    transform(gLag_, gLag_+nx_, mu_x_, gLag_, plus<double>());

    // Number of SQP iterations
    int iter = 0;

    // Number of line-search iterations
    int ls_iter = 0;

    // Last linesearch successfull
    bool ls_success = true;

    // Reset
    merit_mem_.clear();
    sigma_ = 0.;    // NOTE: Move this into the main optimization loop

    // Default stepsize
    double t = 0;

    // MAIN OPTIMIZATION LOOP
    while (true) {

      // Primal infeasability
      double pr_inf = primalInfeasibility(x_, getPtr(lbx), getPtr(ubx),
                                          gk_, getPtr(lbg), getPtr(ubg));

      // inf-norm of lagrange gradient
      double gLag_norminf = casadi_norm_inf(nx_, gLag_);

      // inf-norm of step
      double dx_norminf = casadi_norm_inf(nx_, dx_);

      // Print header occasionally
      if (iter % 10 == 0) printIteration(userOut());

      // Printing information about the actual iterate
      printIteration(userOut(), iter, fk_, pr_inf, gLag_norminf, dx_norminf,
                     reg_, ls_iter, ls_success);

      if (gather_stats_) {
        Dict iterations = stats_["iterations"];
        std::vector<double> tmp=iterations["inf_pr"];
        tmp.push_back(pr_inf);
        iterations["inf_pr"] = tmp;

        tmp=iterations["inf_du"];
        tmp.push_back(gLag_norminf);
        iterations["inf_du"] = tmp;

        tmp=iterations["d_norm"];
        tmp.push_back(dx_norminf);
        iterations["d_norm"] = tmp;

        std::vector<int> tmp2=iterations["ls_trials"];
        tmp2.push_back(ls_iter);
        iterations["ls_trials"] = tmp2;

        tmp=iterations["obj"];
        tmp.push_back(fk_);
        iterations["obj"] = tmp;

        stats_["iterations"] = iterations;
      }

      // Call callback function if present
      if (!fcallback_.isNull()) {
        double time1 = clock();

        if (!output(NLPSOL_F).is_empty()) output(NLPSOL_F).set(fk_);
        if (!output(NLPSOL_X).is_empty()) output(NLPSOL_X).setNZ(x_);
        if (!output(NLPSOL_LAM_G).is_empty()) output(NLPSOL_LAM_G).setNZ(mu_);
        if (!output(NLPSOL_LAM_X).is_empty()) output(NLPSOL_LAM_X).setNZ(mu_x_);
        if (!output(NLPSOL_G).is_empty()) output(NLPSOL_G).setNZ(gk_);

        Dict iteration;
        iteration["iter"] = iter;
        iteration["inf_pr"] = pr_inf;
        iteration["inf_du"] = gLag_norminf;
        iteration["d_norm"] = dx_norminf;
        iteration["ls_trials"] = ls_iter;
        iteration["obj"] = fk_;
        stats_["iteration"] = iteration;

        double time2 = clock();
        t_callback_prepare_ += (time2-time1)/CLOCKS_PER_SEC;
        time1 = clock();

        for (int i=0; i<NLPSOL_NUM_OUT; ++i) {
          fcallback_.setInput(output(i), i);
        }
        fcallback_.evaluate();
        double ret_double;
        fcallback_.getOutput(ret_double);
        int ret = static_cast<int>(ret_double);
        time2 = clock();
        t_callback_fun_ += (time2-time1)/CLOCKS_PER_SEC;
        if (ret) {
          userOut() << endl;
          userOut() << "casadi::SQPMethod: aborted by callback..." << endl;
          stats_["return_status"] = "User_Requested_Stop";
          break;
        }
      }

      // Checking convergence criteria
      if (pr_inf < tol_pr_ && gLag_norminf < tol_du_) {
        userOut() << endl;
        userOut() << "casadi::SQPMethod: Convergence achieved after "
                  << iter << " iterations." << endl;
        stats_["return_status"] = "Solve_Succeeded";
        break;
      }

      if (iter >= max_iter_) {
        userOut() << endl;
        userOut() << "casadi::SQPMethod: Maximum number of iterations reached." << endl;
        stats_["return_status"] = "Maximum_Iterations_Exceeded";
        break;
      }

      if (iter > 0 && dx_norminf <= min_step_size_) {
        userOut() << endl;
        userOut() << "casadi::SQPMethod: Search direction becomes too small without "
            "convergence criteria being met." << endl;
        stats_["return_status"] = "Search_Direction_Becomes_Too_Small";
        break;
      }

      // Start a new iteration
      iter++;

      log("Formulating QP");
      // Formulate the QP
      transform(lbx.begin(), lbx.end(), x_, qp_LBX_, minus<double>());
      transform(ubx.begin(), ubx.end(), x_, qp_UBX_, minus<double>());
      transform(lbg.begin(), lbg.end(), gk_, qp_LBA_, minus<double>());
      transform(ubg.begin(), ubg.end(), gk_, qp_UBA_, minus<double>());

      // Solve the QP
      solve_QP(Bk_, gf_, qp_LBX_, qp_UBX_, Jk_, qp_LBA_,
               qp_UBA_, dx_, qp_DUAL_X_, qp_DUAL_A_);
      log("QP solved");

      // Detecting indefiniteness
      double gain = casadi_quad_form(Bk_, Hsp_, dx_);
      if (gain < 0) {
        casadi_warning("Indefinite Hessian detected...");
      }

      // Calculate penalty parameter of merit function
      sigma_ = std::max(sigma_, 1.01*casadi_norm_inf(nx_, qp_DUAL_X_));
      sigma_ = std::max(sigma_, 1.01*casadi_norm_inf(ng_, qp_DUAL_A_));

      // Calculate L1-merit function in the actual iterate
      double l1_infeas = primalInfeasibility(x_, getPtr(lbx), getPtr(ubx),
                                             gk_, getPtr(lbg), getPtr(ubg));

      // Right-hand side of Armijo condition
      double F_sens = casadi_inner_prod(nx_, dx_, gf_);
      double L1dir = F_sens - sigma_ * l1_infeas;
      double L1merit = fk_ + sigma_ * l1_infeas;

      // Storing the actual merit function value in a list
      merit_mem_.push_back(L1merit);
      if (merit_mem_.size() > merit_memsize_) {
        merit_mem_.pop_front();
      }
      // Stepsize
      t = 1.0;
      double fk_cand;
      // Merit function value in candidate
      double L1merit_cand = 0;

      // Reset line-search counter, success marker
      ls_iter = 0;
      ls_success = true;

      // Line-search
      log("Starting line-search");
      if (max_iter_ls_>0) { // max_iter_ls_== 0 disables line-search

        // Line-search loop
        while (true) {
          for (int i=0; i<nx_; ++i) x_cand_[i] = x_[i] + t * dx_[i];

          try {
            // Evaluating objective and constraints
            fk_cand = eval_f(x_cand_);
            eval_g(x_cand_, gk_cand_);
          } catch(const CasadiException& ex) {
            // Silent ignore; line-search failed
            ls_iter++;
            // Backtracking
            t = beta_ * t;
            continue;
          }

          ls_iter++;

          // Calculating merit-function in candidate
          l1_infeas = primalInfeasibility(x_cand_, getPtr(lbx), getPtr(ubx),
                                          gk_cand_, getPtr(lbg), getPtr(ubg));

          L1merit_cand = fk_cand + sigma_ * l1_infeas;
          // Calculating maximal merit function value so far
          double meritmax = *max_element(merit_mem_.begin(), merit_mem_.end());
          if (L1merit_cand <= meritmax + t * c1_ * L1dir) {
            // Accepting candidate
            log("Line-search completed, candidate accepted");
            break;
          }

          // Line-search not successful, but we accept it.
          if (ls_iter == max_iter_ls_) {
            ls_success = false;
            log("Line-search completed, maximum number of iterations");
            break;
          }

          // Backtracking
          t = beta_ * t;
        }

        // Candidate accepted, update dual variables
        for (int i=0; i<ng_; ++i) mu_[i] = t * qp_DUAL_A_[i] + (1 - t) * mu_[i];
        for (int i=0; i<nx_; ++i) mu_x_[i] = t * qp_DUAL_X_[i] + (1 - t) * mu_x_[i];

        // Candidate accepted, update the primal variable
        copy(x_, x_+nx_, x_old_);
        copy(x_cand_, x_cand_+nx_, x_);

      } else {
        // Full step
        copy_n(qp_DUAL_A_, ng_, mu_);
        copy_n(qp_DUAL_X_, nx_, mu_x_);

        copy(x_, x_+nx_, x_old_);
        // x+=dx
        transform(x_, x_+nx_, dx_, x_, plus<double>());
      }

      if (!exact_hessian_) {
        // Evaluate the gradient of the Lagrangian with the old x but new mu (for BFGS)
        copy(gf_, gf_+nx_, gLag_old_);
        if (ng_>0) casadi_mv_t(Jk_, Asp_, mu_, gLag_old_);
        // gLag_old += mu_x_;
        transform(gLag_old_, gLag_old_+nx_, mu_x_, gLag_old_, plus<double>());
      }

      // Evaluate the constraint Jacobian
      log("Evaluating jac_g");
      eval_jac_g(x_, gk_, Jk_);

      // Evaluate the gradient of the objective function
      log("Evaluating grad_f");
      eval_grad_f(x_, &fk_, gf_);

      // Evaluate the gradient of the Lagrangian with the new x and new mu
      copy(gf_, gf_+nx_, gLag_);
      if (ng_>0) casadi_mv_t(Jk_, Asp_, mu_, gLag_);

      // gLag += mu_x_;
      transform(gLag_, gLag_+nx_, mu_x_, gLag_, plus<double>());

      // Updating Lagrange Hessian
      if (!exact_hessian_) {
        log("Updating Hessian (BFGS)");
        // BFGS with careful updates and restarts
        if (iter % lbfgs_memory_ == 0) {
          // Reset Hessian approximation by dropping all off-diagonal entries
          const int* colind = Hsp_.colind();      // Access sparsity (column offset)
          int ncol = Hsp_.size2();
          const int* row = Hsp_.row();            // Access sparsity (row)
          for (int cc=0; cc<ncol; ++cc) {     // Loop over the columns of the Hessian
            for (int el=colind[cc]; el<colind[cc+1]; ++el) {
              // Loop over the nonzero elements of the column
              if (cc!=row[el]) Bk_[el] = 0;               // Remove if off-diagonal entries
            }
          }
        }

        // Pass to BFGS update function
        bfgs_.setInputNZ(Bk_, BFGS_BK);
        bfgs_.setInputNZ(x_, BFGS_X);
        bfgs_.setInputNZ(x_old_, BFGS_X_OLD);
        bfgs_.setInputNZ(gLag_, BFGS_GLAG);
        bfgs_.setInputNZ(gLag_old_, BFGS_GLAG_OLD);

        // Update the Hessian approximation
        bfgs_.evaluate();

        // Get the updated Hessian
        bfgs_.getOutputNZ(Bk_);

      } else {
        // Exact Hessian
        log("Evaluating hessian");
        eval_h(x_, mu_, 1.0, Bk_);
      }
    }

    double time2 = clock();
    t_mainloop_ = (time2-time1)/CLOCKS_PER_SEC;

    // Save results to outputs
    output(NLPSOL_F).set(fk_);
    output(NLPSOL_X).setNZ(x_);
    output(NLPSOL_LAM_G).setNZ(mu_);
    output(NLPSOL_LAM_X).setNZ(mu_x_);
    output(NLPSOL_G).setNZ(gk_);

    if (hasOption("print_time") && static_cast<bool>(option("print_time"))) {
      // Write timings
      userOut() << "time spent in eval_f: " << t_eval_f_ << " s.";
      if (n_eval_f_>0)
        userOut() << " (" << n_eval_f_ << " calls, " << (t_eval_f_/n_eval_f_)*1000
                  << " ms. average)";
      userOut() << endl;
      userOut() << "time spent in eval_grad_f: " << t_eval_grad_f_ << " s.";
      if (n_eval_grad_f_>0)
        userOut() << " (" << n_eval_grad_f_ << " calls, "
             << (t_eval_grad_f_/n_eval_grad_f_)*1000 << " ms. average)";
      userOut() << endl;
      userOut() << "time spent in eval_g: " << t_eval_g_ << " s.";
      if (n_eval_g_>0)
        userOut() << " (" << n_eval_g_ << " calls, " << (t_eval_g_/n_eval_g_)*1000
                  << " ms. average)";
      userOut() << endl;
      userOut() << "time spent in eval_jac_g: " << t_eval_jac_g_ << " s.";
      if (n_eval_jac_g_>0)
        userOut() << " (" << n_eval_jac_g_ << " calls, "
             << (t_eval_jac_g_/n_eval_jac_g_)*1000 << " ms. average)";
      userOut() << endl;
      userOut() << "time spent in eval_h: " << t_eval_h_ << " s.";
      if (n_eval_h_>1)
        userOut() << " (" << n_eval_h_ << " calls, " << (t_eval_h_/n_eval_h_)*1000
                  << " ms. average)";
      userOut() << endl;
      userOut() << "time spent in main loop: " << t_mainloop_ << " s." << endl;
      userOut() << "time spent in callback function: " << t_callback_fun_ << " s." << endl;
      userOut() << "time spent in callback preparation: " << t_callback_prepare_ << " s." << endl;
    }

    // Save statistics
    stats_["iter_count"] = iter;

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

    // Get the outputs
    for (int i=0; i<n_out(); ++i) {
      if (res[i] != 0) getOutputNZ(res[i], i);
    }
  }

  void Sqpmethod::printIteration(std::ostream &stream) {
    stream << setw(4)  << "iter";
    stream << setw(15) << "objective";
    stream << setw(10) << "inf_pr";
    stream << setw(10) << "inf_du";
    stream << setw(10) << "||d||";
    stream << setw(7) << "lg(rg)";
    stream << setw(3) << "ls";
    stream << ' ';
    stream << endl;
  }

  void Sqpmethod::printIteration(std::ostream &stream, int iter, double obj,
                                   double pr_inf, double du_inf,
                                   double dx_norm, double rg, int ls_trials, bool ls_success) {
    stream << setw(4) << iter;
    stream << scientific;
    stream << setw(15) << setprecision(6) << obj;
    stream << setw(10) << setprecision(2) << pr_inf;
    stream << setw(10) << setprecision(2) << du_inf;
    stream << setw(10) << setprecision(2) << dx_norm;
    stream << fixed;
    if (rg>0) {
      stream << setw(7) << setprecision(2) << log10(rg);
    } else {
      stream << setw(7) << "-";
    }
    stream << setw(3) << ls_trials;
    stream << (ls_success ? ' ' : 'F');
    stream << endl;
  }

  void Sqpmethod::reset_h() {
    // Initial Hessian approximation of BFGS
    if (!exact_hessian_) {
      copy_n(B_init_.ptr(), Hsp_.nnz(), Bk_);
    }
  }

  double Sqpmethod::getRegularization(const double* H) {
    const int* colind = Hsp_.colind();
    int ncol = Hsp_.size2();
    const int* row = Hsp_.row();
    double reg_param = 0;
    for (int cc=0; cc<ncol; ++cc) {
      double mineig = 0;
      for (int el=colind[cc]; el<colind[cc+1]; ++el) {
        int rr = row[el];
        if (rr == cc) {
          mineig += H[el];
        } else {
          mineig -= fabs(H[el]);
        }
      }
      reg_param = fmin(reg_param, mineig);
    }
    return -reg_param;
  }

  void Sqpmethod::regularize(double* H, double reg) {
    const int* colind = Hsp_.colind();
    int ncol = Hsp_.size2();
    const int* row = Hsp_.row();

    for (int cc=0; cc<ncol; ++cc) {
      for (int el=colind[cc]; el<colind[cc+1]; ++el) {
        int rr = row[el];
        if (rr==cc) {
          H[el] += reg;
        }
      }
    }
  }


  void Sqpmethod::eval_h(const double* x, const double* lambda, double sigma, double* H) {
    try {
      // Get function
      Function& hessLag = this->hessLag();

      // Pass the argument to the function
      hessLag.setInputNZ(x, HESSLAG_X);
      hessLag.setInputNZ(p_, HESSLAG_P);
      hessLag.setInputNZ(&sigma, HESSLAG_LAM_F);
      hessLag.setInputNZ(lambda, HESSLAG_LAM_G);

      // Evaluate
      hessLag.evaluate();

      // Get results
      hessLag.getOutputNZ(H);

      // Determing regularization parameter with Gershgorin theorem
      if (regularize_) {
        reg_ = getRegularization(H);
        if (reg_ > 0) {
          regularize(H, reg_);
        }
      }

    } catch(exception& ex) {
      userOut<true, PL_WARN>() << "eval_h failed: " << ex.what() << endl;
      throw;
    }
  }

  void Sqpmethod::eval_g(const double* x, double* g) {
    try {
      double time1 = clock();

      // Quick return if no constraints
      if (ng_==0) return;

      // Pass the argument to the function
      nlp_.setInputNZ(x, NL_X);
      nlp_.setInputNZ(p_, NL_P);

      // Evaluate the function and tape
      nlp_.evaluate();

      // Ge the result
      nlp_.getOutputNZ(g, NL_G);

      double time2 = clock();
      t_eval_g_ += (time2-time1)/CLOCKS_PER_SEC;
      n_eval_g_ += 1;
    } catch(exception& ex) {
      userOut<true, PL_WARN>() << "eval_g failed: " << ex.what() << endl;
      throw;
    }
  }

  void Sqpmethod::eval_jac_g(const double* x, double* g, double* J) {
    try {
      double time1 = clock();

      // Quich finish if no constraints
      if (ng_==0) return;

      // Get function
      Function& jacG = this->jacG();

      // Pass the argument to the function
      jacG.setInputNZ(x, NL_X);
      jacG.setInputNZ(p_, NL_P);

      // Evaluate the function
      jacG.evaluate();

      // Get the output
      jacG.getOutputNZ(g, 1+NL_G);
      jacG.getOutputNZ(J);

      double time2 = clock();
      t_eval_jac_g_ += (time2-time1)/CLOCKS_PER_SEC;
      n_eval_jac_g_ += 1;

    } catch(exception& ex) {
      userOut<true, PL_WARN>() << "eval_jac_g failed: " << ex.what() << endl;
      throw;
    }
  }

  void Sqpmethod::eval_grad_f(const double* x, double* f, double* grad_f) {
    try {
      double time1 = clock();

      // Get function
      Function& gradF = this->gradF();

      // Pass the argument to the function
      gradF.setInputNZ(x, NL_X);
      gradF.setInputNZ(p_, NL_P);

      // Evaluate, adjoint mode
      gradF.evaluate();

      // Get the result
      gradF.getOutputNZ(grad_f);
      gradF.getOutputNZ(f, 1+NL_X);

      double time2 = clock();
      t_eval_grad_f_ += (time2-time1)/CLOCKS_PER_SEC;
      n_eval_grad_f_ += 1;

    } catch(exception& ex) {
      userOut<true, PL_WARN>() << "eval_grad_f failed: " << ex.what() << endl;
      throw;
    }
  }

  double Sqpmethod::eval_f(const double* x) {
    try {
       // Log time
      double time1 = clock();

      // Pass the argument to the function
      nlp_.setInputNZ(x, NL_X);
      nlp_.setInputNZ(p_, NL_P);

      // Evaluate the function
      nlp_.evaluate();

      // Get the result
      double f;
      nlp_.getOutputNZ(&f, NL_F);

      double time2 = clock();
      t_eval_f_ += (time2-time1)/CLOCKS_PER_SEC;
      n_eval_f_ += 1;

      return f;
    } catch(exception& ex) {
      userOut<true, PL_WARN>() << "eval_f failed: " << ex.what() << endl;
      throw;
    }
  }

  void Sqpmethod::solve_QP(const double* H, const double* g, const double* lbx, const double* ubx,
                           const double* A, const double* lbA, const double* ubA,
                           double* x_opt, double* lambda_x_opt, double* lambda_A_opt) {

    // Pass data to QP solver
    qpsol_.setInputNZ(H, QPSOL_H);
    qpsol_.setInputNZ(g, QPSOL_G);

    // Hot-starting if possible
    qpsol_.setInputNZ(x_opt, QPSOL_X0);

    // Pass simple bounds
    qpsol_.setInputNZ(lbx, QPSOL_LBX);
    qpsol_.setInputNZ(ubx, QPSOL_UBX);

    // Pass linear bounds
    if (ng_>0) {
      qpsol_.setInputNZ(A, QPSOL_A);
      qpsol_.setInputNZ(lbA, QPSOL_LBA);
      qpsol_.setInputNZ(ubA, QPSOL_UBA);
    }

    // Solve the QP
    qpsol_.evaluate();

    // Get the optimal solution
    qpsol_.getOutputNZ(x_opt, QPSOL_X);
    qpsol_.getOutputNZ(lambda_x_opt, QPSOL_LAM_X);
    qpsol_.getOutputNZ(lambda_A_opt, QPSOL_LAM_A);
  }

  double Sqpmethod::primalInfeasibility(const double* x, const double* lbx, const double* ubx,
                                        const double* g, const double* lbg, const double* ubg) {
    // Linf-norm of the primal infeasibility
    double pr_inf = 0;

    // Bound constraints
    for (int j=0; j<nx_; ++j) {
      pr_inf = fmax(pr_inf, lbx[j] - x[j]);
      pr_inf = fmax(pr_inf, x[j] - ubx[j]);
    }

    // Nonlinear constraints
    for (int j=0; j<ng_; ++j) {
      pr_inf = fmax(pr_inf, lbg[j] - g[j]);
      pr_inf = fmax(pr_inf, g[j] - ubg[j]);
    }

    return pr_inf;
  }

} // namespace casadi
