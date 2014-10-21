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
  int CASADI_NLPSOLVER_SQPMETHOD_EXPORT
      casadi_register_nlpsolver_sqpmethod(NlpSolverInternal::Plugin* plugin) {
    plugin->creator = Sqpmethod::creator;
    plugin->name = "sqpmethod";
    plugin->doc = Sqpmethod::meta_doc.c_str();
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOLVER_SQPMETHOD_EXPORT casadi_load_nlpsolver_sqpmethod() {
    NlpSolverInternal::registerPlugin(casadi_register_nlpsolver_sqpmethod);
  }

  Sqpmethod::Sqpmethod(const Function& nlp) : NlpSolverInternal(nlp) {
    casadi_warning("The SQP method is under development");
    addOption("qp_solver",         OT_STRING,   GenericType(),
              "The QP solver to be used by the SQP method");
    addOption("qp_solver_options", OT_DICTIONARY, GenericType(),
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

    std::string qp_solver_name = getOption("qp_solver");
    qp_solver_ = QpSolver(qp_solver_name,
                          qpStruct("h", H_sparsity, "a", A_sparsity));

    // Set options if provided
    if (hasSetOption("qp_solver_options")) {
      Dictionary qp_solver_options = getOption("qp_solver_options");
      qp_solver_.setOption(qp_solver_options);
    }
    qp_solver_.init();

    // Lagrange multipliers of the NLP
    mu_.resize(ng_);
    mu_x_.resize(nx_);

    // Lagrange gradient in the next iterate
    gLag_.resize(nx_);
    gLag_old_.resize(nx_);

    // Current linearization point
    x_.resize(nx_);
    x_cand_.resize(nx_);
    x_old_.resize(nx_);

    // Constraint function value
    gk_.resize(ng_);
    gk_cand_.resize(ng_);

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

    // Gradient of the objective
    gf_.resize(nx_);

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
      cout << "This is casadi::SQPMethod." << endl;
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

  void Sqpmethod::evaluate() {
    if (inputs_check_) checkInputs();
    checkInitialBounds();

    if (gather_stats_) {
      Dictionary iterations;
      iterations["inf_pr"] = std::vector<double>();
      iterations["inf_du"] = std::vector<double>();
      iterations["ls_trials"] = std::vector<double>();
      iterations["d_norm"] = std::vector<double>();
      iterations["obj"] = std::vector<double>();
      stats_["iterations"] = iterations;
    }


    // Get problem data
    const vector<double>& x_init = input(NLP_SOLVER_X0).data();
    const vector<double>& lbx = input(NLP_SOLVER_LBX).data();
    const vector<double>& ubx = input(NLP_SOLVER_UBX).data();
    const vector<double>& lbg = input(NLP_SOLVER_LBG).data();
    const vector<double>& ubg = input(NLP_SOLVER_UBG).data();

    // Set linearization point to initial guess
    copy(x_init.begin(), x_init.end(), x_.begin());

    // Initialize Lagrange multipliers of the NLP
    copy(input(NLP_SOLVER_LAM_G0).begin(), input(NLP_SOLVER_LAM_G0).end(), mu_.begin());
    copy(input(NLP_SOLVER_LAM_X0).begin(), input(NLP_SOLVER_LAM_X0).end(), mu_x_.begin());

    t_eval_f_ = t_eval_grad_f_ = t_eval_g_ = t_eval_jac_g_ = t_eval_h_ =
        t_callback_fun_ = t_callback_prepare_ = t_mainloop_ = 0;

    n_eval_f_ = n_eval_grad_f_ = n_eval_g_ = n_eval_jac_g_ = n_eval_h_ = 0;

    double time1 = clock();

    // Initial constraint Jacobian
    eval_jac_g(x_, gk_, Jk_);

    // Initial objective gradient
    eval_grad_f(x_, fk_, gf_);

    // Initialize or reset the Hessian or Hessian approximation
    reg_ = 0;
    if (exact_hessian_) {
      eval_h(x_, mu_, 1.0, Bk_);
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
    sigma_ = 0.;    // NOTE: Move this into the main optimization loop

    // Default stepsize
    double t = 0;

    // MAIN OPTIMIZATION LOOP
    while (true) {

      // Primal infeasability
      double pr_inf = primalInfeasibility(x_, lbx, ubx, gk_, lbg, ubg);

      // inf-norm of lagrange gradient
      double gLag_norminf = norm_inf(gLag_);

      // inf-norm of step
      double dx_norminf = norm_inf(dx_);

      // Print header occasionally
      if (iter % 10 == 0) printIteration(cout);

      // Printing information about the actual iterate
      printIteration(cout, iter, fk_, pr_inf, gLag_norminf, dx_norminf, reg_, ls_iter, ls_success);

      if (gather_stats_) {
        Dictionary & iterations = stats_["iterations"];
        static_cast<std::vector<double> &>(iterations["inf_pr"]).push_back(pr_inf);
        static_cast<std::vector<double> &>(iterations["inf_du"]).push_back(gLag_norminf);
        static_cast<std::vector<double> &>(iterations["d_norm"]).push_back(dx_norminf);
        static_cast<std::vector<double> &>(iterations["ls_trials"]).push_back(ls_iter);
        static_cast<std::vector<double> &>(iterations["obj"]).push_back(fk_);
      }

      // Call callback function if present
      if (!callback_.isNull()) {
        double time1 = clock();

        if (!output(NLP_SOLVER_F).isEmpty()) output(NLP_SOLVER_F).set(fk_);
        if (!output(NLP_SOLVER_X).isEmpty()) output(NLP_SOLVER_X).set(x_);
        if (!output(NLP_SOLVER_LAM_G).isEmpty()) output(NLP_SOLVER_LAM_G).set(mu_);
        if (!output(NLP_SOLVER_LAM_X).isEmpty()) output(NLP_SOLVER_LAM_X).set(mu_x_);
        if (!output(NLP_SOLVER_G).isEmpty()) output(NLP_SOLVER_G).set(gk_);

        Dictionary iteration;
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
        int ret = callback_(ref_, user_data_);
        time2 = clock();
        t_callback_fun_ += (time2-time1)/CLOCKS_PER_SEC;
        if (ret) {
          cout << endl;
          cout << "casadi::SQPMethod: aborted by callback..." << endl;
          stats_["return_status"] = "User_Requested_Stop";
          break;
        }
      }

      // Checking convergence criteria
      if (pr_inf < tol_pr_ && gLag_norminf < tol_du_) {
        cout << endl;
        cout << "casadi::SQPMethod: Convergence achieved after " << iter << " iterations." << endl;
        stats_["return_status"] = "Solve_Succeeded";
        break;
      }

      if (iter >= max_iter_) {
        cout << endl;
        cout << "casadi::SQPMethod: Maximum number of iterations reached." << endl;
        stats_["return_status"] = "Maximum_Iterations_Exceeded";
        break;
      }

      if (iter > 0 && dx_norminf <= min_step_size_) {
        cout << endl;
        cout << "casadi::SQPMethod: Search direction becomes too small without "
            "convergence criteria being met." << endl;
        stats_["return_status"] = "Search_Direction_Becomes_Too_Small";
        break;
      }

      // Start a new iteration
      iter++;

      log("Formulating QP");
      // Formulate the QP
      transform(lbx.begin(), lbx.end(), x_.begin(), qp_LBX_.begin(), minus<double>());
      transform(ubx.begin(), ubx.end(), x_.begin(), qp_UBX_.begin(), minus<double>());
      transform(lbg.begin(), lbg.end(), gk_.begin(), qp_LBA_.begin(), minus<double>());
      transform(ubg.begin(), ubg.end(), gk_.begin(), qp_UBA_.begin(), minus<double>());

      // Solve the QP
      solve_QP(Bk_, gf_, qp_LBX_, qp_UBX_, Jk_, qp_LBA_, qp_UBA_, dx_, qp_DUAL_X_, qp_DUAL_A_);
      log("QP solved");

      // Detecting indefiniteness
      double gain = quad_form(dx_, Bk_);
      if (gain < 0) {
        casadi_warning("Indefinite Hessian detected...");
      }

      // Calculate penalty parameter of merit function
      sigma_ = std::max(sigma_, 1.01*norm_inf(qp_DUAL_X_));
      sigma_ = std::max(sigma_, 1.01*norm_inf(qp_DUAL_A_));

      // Calculate L1-merit function in the actual iterate
      double l1_infeas = primalInfeasibility(x_, lbx, ubx, gk_, lbg, ubg);

      // Right-hand side of Armijo condition
      double F_sens = inner_prod(dx_, gf_);
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
            eval_f(x_cand_, fk_cand);
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
          l1_infeas = primalInfeasibility(x_cand_, lbx, ubx, gk_cand_, lbg, ubg);

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
        copy(x_.begin(), x_.end(), x_old_.begin());
        copy(x_cand_.begin(), x_cand_.end(), x_.begin());

      } else {
        // Full step
        copy(qp_DUAL_A_.begin(), qp_DUAL_A_.end(), mu_.begin());
        copy(qp_DUAL_X_.begin(), qp_DUAL_X_.end(), mu_x_.begin());

        copy(x_.begin(), x_.end(), x_old_.begin());
        // x+=dx
        transform(x_.begin(), x_.end(), dx_.begin(), x_.begin(), plus<double>());
      }

      if (!exact_hessian_) {
        // Evaluate the gradient of the Lagrangian with the old x but new mu (for BFGS)
        copy(gf_.begin(), gf_.end(), gLag_old_.begin());
        if (ng_>0) DMatrix::mul_no_alloc_tn(Jk_, mu_, gLag_old_);
        // gLag_old += mu_x_;
        transform(gLag_old_.begin(), gLag_old_.end(), mu_x_.begin(), gLag_old_.begin(),
                  plus<double>());
      }

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
          for (int cc=0; cc<colind.size()-1; ++cc) {     // Loop over the columns of the Hessian
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
        if (monitored("bfgs")) {
          cout << "x = " << x_ << endl;
          cout << "BFGS = "  << endl;
          Bk_.printSparse();
        }
      } else {
        // Exact Hessian
        log("Evaluating hessian");
        eval_h(x_, mu_, 1.0, Bk_);
      }
    }

    double time2 = clock();
    t_mainloop_ = (time2-time1)/CLOCKS_PER_SEC;

    // Save results to outputs
    output(NLP_SOLVER_F).set(fk_);
    output(NLP_SOLVER_X).set(x_);
    output(NLP_SOLVER_LAM_G).set(mu_);
    output(NLP_SOLVER_LAM_X).set(mu_x_);
    output(NLP_SOLVER_G).set(gk_);

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

  double Sqpmethod::quad_form(const std::vector<double>& x, const DMatrix& A) {
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
        // Get row
        int rr = A_row[el];

        // Add contribution
        ret += x[cc]*A_data[el]*x[rr];
      }
    }

    return ret;
  }

  void Sqpmethod::reset_h() {
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

  double Sqpmethod::getRegularization(const Matrix<double>& H) {
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

  void Sqpmethod::regularize(Matrix<double>& H, double reg) {
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


  void Sqpmethod::eval_h(const std::vector<double>& x, const std::vector<double>& lambda,
                           double sigma, Matrix<double>& H) {
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

  void Sqpmethod::eval_g(const std::vector<double>& x, std::vector<double>& g) {
    try {
      double time1 = clock();

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

      double time2 = clock();
      t_eval_g_ += (time2-time1)/CLOCKS_PER_SEC;
      n_eval_g_ += 1;
    } catch(exception& ex) {
      cerr << "eval_g failed: " << ex.what() << endl;
      throw;
    }
  }

  void Sqpmethod::eval_jac_g(const std::vector<double>& x, std::vector<double>& g,
                               Matrix<double>& J) {
    try {
      double time1 = clock();

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

      double time2 = clock();
      t_eval_jac_g_ += (time2-time1)/CLOCKS_PER_SEC;
      n_eval_jac_g_ += 1;

    } catch(exception& ex) {
      cerr << "eval_jac_g failed: " << ex.what() << endl;
      throw;
    }
  }

  void Sqpmethod::eval_grad_f(const std::vector<double>& x, double& f,
                                std::vector<double>& grad_f) {
    try {
      double time1 = clock();

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
      double time2 = clock();
      t_eval_grad_f_ += (time2-time1)/CLOCKS_PER_SEC;
      n_eval_grad_f_ += 1;

    } catch(exception& ex) {
      cerr << "eval_grad_f failed: " << ex.what() << endl;
      throw;
    }
  }

  void Sqpmethod::eval_f(const std::vector<double>& x, double& f) {
    try {
       // Log time
      double time1 = clock();

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
      double time2 = clock();
      t_eval_f_ += (time2-time1)/CLOCKS_PER_SEC;
      n_eval_f_ += 1;

    } catch(exception& ex) {
      cerr << "eval_f failed: " << ex.what() << endl;
      throw;
    }
  }

  void Sqpmethod::solve_QP(const Matrix<double>& H, const std::vector<double>& g,
                             const std::vector<double>& lbx, const std::vector<double>& ubx,
                             const Matrix<double>& A, const std::vector<double>& lbA,
                             const std::vector<double>& ubA,
                             std::vector<double>& x_opt, std::vector<double>& lambda_x_opt,
                             std::vector<double>& lambda_A_opt) {

    // Pass data to QP solver
    qp_solver_.setInput(H, QP_SOLVER_H);
    qp_solver_.setInput(g, QP_SOLVER_G);

    // Hot-starting if possible
    qp_solver_.setInput(x_opt, QP_SOLVER_X0);

    //TODO(Joel): Fix hot-starting of dual variables
    //qp_solver_.setInput(lambda_A_opt, QP_SOLVER_LAMBDA_INIT);

    // Pass simple bounds
    qp_solver_.setInput(lbx, QP_SOLVER_LBX);
    qp_solver_.setInput(ubx, QP_SOLVER_UBX);

    // Pass linear bounds
    if (ng_>0) {
      qp_solver_.setInput(A, QP_SOLVER_A);
      qp_solver_.setInput(lbA, QP_SOLVER_LBA);
      qp_solver_.setInput(ubA, QP_SOLVER_UBA);
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
    qp_solver_.evaluate();

    // Get the optimal solution
    qp_solver_.getOutput(x_opt, QP_SOLVER_X);
    qp_solver_.getOutput(lambda_x_opt, QP_SOLVER_LAM_X);
    qp_solver_.getOutput(lambda_A_opt, QP_SOLVER_LAM_A);
    if (monitored("dx")) {
      cout << "dx = " << x_opt << endl;
    }
  }

  double Sqpmethod::primalInfeasibility(const std::vector<double>& x,
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

} // namespace casadi
