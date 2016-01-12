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
#include "casadi/core/calculus.hpp"

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
    addOption("max_iter",           OT_INT,      50,
              "Maximum number of SQP iterations");
    addOption("max_iter_ls",        OT_INT,       3,
              "Maximum number of linesearch iterations");
    addOption("tol_pr",            OT_DOUBLE,       1e-6,
              "Stopping criterion for primal infeasibility");
    addOption("tol_du",            OT_DOUBLE,       1e-6,
              "Stopping criterion for dual infeasability");
    addOption("c1",                OT_DOUBLE,       1E-4,
              "Armijo condition, coefficient of decrease in merit");
    addOption("beta",              OT_DOUBLE,       0.8,
              "Line-search parameter, restoration factor of stepsize");
    addOption("merit_memory",      OT_INT,      4,
              "Size of memory to store history of merit function values");
    addOption("lbfgs_memory",      OT_INT,     10,
              "Size of L-BFGS memory.");
    addOption("regularize",        OT_BOOL,  false,
              "Automatic regularization of Lagrange Hessian.");
    addOption("print_header",      OT_BOOL,   true,
              "Print the header with problem statistics");
    addOption("min_step_size",     OT_DOUBLE,   1e-10,
              "The size (inf-norm) of the step size should not become smaller than this.");

    // Monitors
    addOption("monitor",      OT_STRINGVECTOR, GenericType(),  "",
              "eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h|qp|dx|bfgs", true);
    addOption("print_time",         OT_BOOL,       true,
              "Print information about execution time");
  }


  Sqpmethod::~Sqpmethod() {
  }

  void Sqpmethod::init(const Dict& opts) {
    // Call the init method of the base class
    Nlpsol::init(opts);

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
    setup_f();
    setup_g();
    setup_grad_f();
    setup_jac_g();
    if (exact_hessian_) {
      setup_hess_l(false, true, true);
    }

    // Allocate a QP solver
    Hsp_ = exact_hessian_ ? hesslag_sp_ : Sparsity::dense(nx_, nx_);
    Asp_ = jac_g_fcn_.isNull() ? Sparsity(0, nx_) : jac_g_fcn_.sparsity_out(1);

    // QP solver options
    Dict qpsol_options;
    if (hasSetOption("qpsol_options")) {
      qpsol_options = option("qpsol_options");
    }

    // Allocate a QP solver
    qpsol_ = qpsol("qpsol", option("qpsol"),
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
      SX qk = mtimes(Bk, sk);

      // Calculating theta
      SX skBksk = dot(sk, qk);
      SX omega = if_else(dot(yk, sk) < 0.2 * dot(sk, qk),
                         0.8 * skBksk / (skBksk - dot(sk, yk)),
                         1);
      yk = omega * yk + (1 - omega) * qk;
      SX theta = 1. / dot(sk, yk);
      SX phi = 1. / dot(qk, sk);
      SX Bk_new = Bk + theta * mtimes(yk, yk.T()) - phi * mtimes(qk, qk.T());

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
      B_init_ = project(DM::eye(nx_), Hsp_);
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

  void Sqpmethod::set_work(Memory& mem, const double**& arg, double**& res,
                                int*& iw, double*& w) const {
    SqpmethodMemory& m = dynamic_cast<SqpmethodMemory&>(mem);

    // Set work in base classes
    Nlpsol::set_work(mem, arg, res, iw, w);

    // Lagrange multipliers of the NLP
    m.mu = w; w += ng_;
    m.mu_x = w; w += nx_;

    // Current linearization point
    m.xk = w; w += nx_;
    m.x_cand = w; w += nx_;
    m.x_old = w; w += nx_;

    // Lagrange gradient in the next iterate
    m.gLag = w; w += nx_;
    m.gLag_old = w; w += nx_;

    // Constraint function value
    m.gk = w; w += ng_;
    m.gk_cand = w; w += ng_;

    // Gradient of the objective
    m.gf = w; w += nx_;

    // Bounds of the QP
    m.qp_LBA = w; w += ng_;
    m.qp_UBA = w; w += ng_;
    m.qp_LBX = w; w += nx_;
    m.qp_UBX = w; w += nx_;

    // QP solution
    m.dx = w; w += nx_;
    m.qp_DUAL_X = w; w += nx_;
    m.qp_DUAL_A = w; w += ng_;

    // Hessian approximation
    m.Bk = w; w += Hsp_.nnz();

    // Jacobian
    m.Jk = w; w += Asp_.nnz();
  }

  void Sqpmethod::solve(Memory& mem) const {
    SqpmethodMemory& m = dynamic_cast<SqpmethodMemory&>(mem);

    // Check the provided inputs
    checkInputs(mem);

    // Set linearization point to initial guess
    casadi_copy(m.x0, nx_, m.xk);

    // Initialize Lagrange multipliers of the NLP
    casadi_copy(m.lam_g0, ng_, m.mu);
    casadi_copy(m.lam_x0, nx_, m.mu_x);

    m.t_eval_f = m.t_eval_grad_f = m.t_eval_g = m.t_eval_jac_g = m.t_eval_h =
      m.t_callback_fun = m.t_callback_prepare = m.t_mainloop = 0;

    m.n_eval_f = m.n_eval_grad_f = m.n_eval_g = m.n_eval_jac_g = m.n_eval_h = 0;

    double time1 = clock();

    // Initial constraint Jacobian
    eval_jac_g(m, m.xk, m.gk, m.Jk);

    // Initial objective gradient
    eval_grad_f(m, m.xk, &m.fk, m.gf);

    // Initialize or reset the Hessian or Hessian approximation
    m.reg = 0;
    if (exact_hessian_) {
      eval_h(m, m.xk, m.mu, 1.0, m.Bk);
    } else {
      reset_h(m);
    }

    // Evaluate the initial gradient of the Lagrangian
    casadi_copy(m.gf, nx_, m.gLag);
    if (ng_>0) casadi_mv(m.Jk, Asp_, m.mu, m.gLag, true);
    // gLag += mu_x_;
    transform(m.gLag, m.gLag+nx_, m.mu_x, m.gLag, plus<double>());

    // Number of SQP iterations
    int iter = 0;

    // Number of line-search iterations
    int ls_iter = 0;

    // Last linesearch successfull
    bool ls_success = true;

    // Reset
    m.merit_mem.clear();
    m.sigma = 0.;    // NOTE: Move this into the main optimization loop

    // Default stepsize
    double t = 0;

    // MAIN OPTIMIZATION LOOP
    while (true) {

      // Primal infeasability
      double pr_inf = primalInfeasibility(m.xk, m.lbx, m.ubx, m.gk, m.lbg, m.ubg);

      // inf-norm of lagrange gradient
      double gLag_norminf = casadi_norm_inf(nx_, m.gLag);

      // inf-norm of step
      double dx_norminf = casadi_norm_inf(nx_, m.dx);

      // Print header occasionally
      if (iter % 10 == 0) printIteration(userOut());

      // Printing information about the actual iterate
      printIteration(userOut(), iter, m.fk, pr_inf, gLag_norminf, dx_norminf,
                     m.reg, ls_iter, ls_success);

      // Call callback function if present
      if (!fcallback_.isNull()) {
        double time1 = clock();

        // Callback inputs
        fill_n(m.arg, fcallback_.n_in(), nullptr);
        m.arg[NLPSOL_F] = &m.fk;
        m.arg[NLPSOL_X] = m.x;
        m.arg[NLPSOL_LAM_G] = m.lam_g;
        m.arg[NLPSOL_LAM_X] = m.lam_x;
        m.arg[NLPSOL_G] = m.g;

        // Callback outputs
        fill_n(m.res, fcallback_.n_out(), nullptr);
        double ret;
        m.arg[0] = &ret;

        // Evaluate
        fcallback_(m.arg, m.res, m.iw, m.w, 0);

        double time2 = clock();
        m.t_callback_fun += (time2-time1)/CLOCKS_PER_SEC;
        if (static_cast<int>(ret)) {
          userOut() << endl;
          userOut() << "casadi::SQPMethod: aborted by callback..." << endl;
          m.return_status = "User_Requested_Stop";
          break;
        }
      }

      // Checking convergence criteria
      if (pr_inf < tol_pr_ && gLag_norminf < tol_du_) {
        userOut() << endl;
        userOut() << "casadi::SQPMethod: Convergence achieved after "
                  << iter << " iterations." << endl;
        m.return_status = "Solve_Succeeded";
        break;
      }

      if (iter >= max_iter_) {
        userOut() << endl;
        userOut() << "casadi::SQPMethod: Maximum number of iterations reached." << endl;
        m.return_status = "Maximum_Iterations_Exceeded";
        break;
      }

      if (iter > 0 && dx_norminf <= min_step_size_) {
        userOut() << endl;
        userOut() << "casadi::SQPMethod: Search direction becomes too small without "
            "convergence criteria being met." << endl;
        m.return_status = "Search_Direction_Becomes_Too_Small";
        break;
      }

      // Start a new iteration
      iter++;

      log("Formulating QP");
      // Formulate the QP
      transform(m.lbx, m.lbx + nx_, m.xk, m.qp_LBX, minus<double>());
      transform(m.ubx, m.ubx + nx_, m.xk, m.qp_UBX, minus<double>());
      transform(m.lbg, m.lbg + ng_, m.gk, m.qp_LBA, minus<double>());
      transform(m.ubg, m.ubg + ng_, m.gk, m.qp_UBA, minus<double>());

      // Solve the QP
      solve_QP(m, m.Bk, m.gf, m.qp_LBX, m.qp_UBX, m.Jk, m.qp_LBA,
               m.qp_UBA, m.dx, m.qp_DUAL_X, m.qp_DUAL_A);
      log("QP solved");

      // Detecting indefiniteness
      double gain = casadi_bilin(m.Bk, Hsp_, m.dx, m.dx);
      if (gain < 0) {
        casadi_warning("Indefinite Hessian detected...");
      }

      // Calculate penalty parameter of merit function
      m.sigma = std::max(m.sigma, 1.01*casadi_norm_inf(nx_, m.qp_DUAL_X));
      m.sigma = std::max(m.sigma, 1.01*casadi_norm_inf(ng_, m.qp_DUAL_A));

      // Calculate L1-merit function in the actual iterate
      double l1_infeas = primalInfeasibility(m.xk, m.lbx, m.ubx, m.gk, m.lbg, m.ubg);

      // Right-hand side of Armijo condition
      double F_sens = casadi_dot(nx_, m.dx, m.gf);
      double L1dir = F_sens - m.sigma * l1_infeas;
      double L1merit = m.fk + m.sigma * l1_infeas;

      // Storing the actual merit function value in a list
      m.merit_mem.push_back(L1merit);
      if (m.merit_mem.size() > merit_memsize_) {
        m.merit_mem.pop_front();
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
          for (int i=0; i<nx_; ++i) m.x_cand[i] = m.xk[i] + t * m.dx[i];

          try {
            // Evaluating objective and constraints
            fk_cand = eval_f(m, m.x_cand);
            eval_g(m, m.x_cand, m.gk_cand);
          } catch(const CasadiException& ex) {
            // Silent ignore; line-search failed
            ls_iter++;
            // Backtracking
            t = beta_ * t;
            continue;
          }

          ls_iter++;

          // Calculating merit-function in candidate
          l1_infeas = primalInfeasibility(m.x_cand, m.lbx, m.ubx, m.gk_cand, m.lbg, m.ubg);

          L1merit_cand = fk_cand + m.sigma * l1_infeas;
          // Calculating maximal merit function value so far
          double meritmax = *max_element(m.merit_mem.begin(), m.merit_mem.end());
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
        for (int i=0; i<ng_; ++i) m.mu[i] = t * m.qp_DUAL_A[i] + (1 - t) * m.mu[i];
        for (int i=0; i<nx_; ++i) m.mu_x[i] = t * m.qp_DUAL_X[i] + (1 - t) * m.mu_x[i];

        // Candidate accepted, update the primal variable
        casadi_copy(m.xk, nx_, m.x_old);
        casadi_copy(m.x_cand, nx_, m.xk);

      } else {
        // Full step
        casadi_copy(m.qp_DUAL_A, ng_, m.mu);
        casadi_copy(m.qp_DUAL_X, nx_, m.mu_x);
        casadi_copy(m.xk, nx_, m.x_old);
        // x+=dx
        transform(m.xk, m.xk+nx_, m.dx, m.xk, plus<double>());
      }

      if (!exact_hessian_) {
        // Evaluate the gradient of the Lagrangian with the old x but new mu (for BFGS)
        casadi_copy(m.gf, nx_, m.gLag_old);
        if (ng_>0) casadi_mv(m.Jk, Asp_, m.mu, m.gLag_old, true);
        // gLag_old += mu_x_;
        transform(m.gLag_old, m.gLag_old+nx_, m.mu_x, m.gLag_old, plus<double>());
      }

      // Evaluate the constraint Jacobian
      log("Evaluating jac_g");
      eval_jac_g(m, m.xk, m.gk, m.Jk);

      // Evaluate the gradient of the objective function
      log("Evaluating grad_f");
      eval_grad_f(m, m.xk, &m.fk, m.gf);

      // Evaluate the gradient of the Lagrangian with the new x and new mu
      casadi_copy(m.gf, nx_, m.gLag);
      if (ng_>0) casadi_mv(m.Jk, Asp_, m.mu, m.gLag, true);

      // gLag += mu_x_;
      transform(m.gLag, m.gLag+nx_, m.mu_x, m.gLag, plus<double>());

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
              if (cc!=row[el]) m.Bk[el] = 0;               // Remove if off-diagonal entries
            }
          }
        }

        // Update the Hessian approximation
        fill_n(m.arg, bfgs_.n_in(), nullptr);
        m.arg[BFGS_BK] = m.Bk;
        m.arg[BFGS_X] = m.xk;
        m.arg[BFGS_X_OLD] = m.x_old;
        m.arg[BFGS_GLAG] = m.gLag;
        m.arg[BFGS_GLAG_OLD] = m.gLag_old;
        fill_n(m.res, bfgs_.n_out(), nullptr);
        m.res[0] = m.Bk;
        bfgs_(m.arg, m.res, m.iw, m.w, 0);

      } else {
        // Exact Hessian
        log("Evaluating hessian");
        eval_h(m, m.xk, m.mu, 1.0, m.Bk);
      }
    }

    double time2 = clock();
    m.t_mainloop = (time2-time1)/CLOCKS_PER_SEC;

    // Save results to outputs
    if (m.f) *m.f = m.fk;
    if (m.x) casadi_copy(m.xk, nx_, m.x);
    if (m.lam_g) casadi_copy(m.mu, ng_, m.lam_g);
    if (m.lam_x) casadi_copy(m.mu_x, nx_, m.lam_x);
    if (m.g) casadi_copy(m.gk, ng_, m.g);

    if (hasOption("print_time") && static_cast<bool>(option("print_time"))) {
      // Write timings
      userOut() << "time spent in eval_f: " << m.t_eval_f << " s.";
      if (m.n_eval_f>0)
        userOut() << " (" << m.n_eval_f << " calls, " << (m.t_eval_f/m.n_eval_f)*1000
                  << " ms. average)";
      userOut() << endl;
      userOut() << "time spent in eval_grad_f: " << m.t_eval_grad_f << " s.";
      if (m.n_eval_grad_f>0)
        userOut() << " (" << m.n_eval_grad_f << " calls, "
             << (m.t_eval_grad_f/m.n_eval_grad_f)*1000 << " ms. average)";
      userOut() << endl;
      userOut() << "time spent in eval_g: " << m.t_eval_g << " s.";
      if (m.n_eval_g>0)
        userOut() << " (" << m.n_eval_g << " calls, " << (m.t_eval_g/m.n_eval_g)*1000
                  << " ms. average)";
      userOut() << endl;
      userOut() << "time spent in eval_jac_g: " << m.t_eval_jac_g << " s.";
      if (m.n_eval_jac_g>0)
        userOut() << " (" << m.n_eval_jac_g << " calls, "
             << (m.t_eval_jac_g/m.n_eval_jac_g)*1000 << " ms. average)";
      userOut() << endl;
      userOut() << "time spent in eval_h: " << m.t_eval_h << " s.";
      if (m.n_eval_h>1)
        userOut() << " (" << m.n_eval_h << " calls, " << (m.t_eval_h/m.n_eval_h)*1000
                  << " ms. average)";
      userOut() << endl;
      userOut() << "time spent in main loop: " << m.t_mainloop << " s." << endl;
      userOut() << "time spent in callback function: " << m.t_callback_fun << " s." << endl;
      userOut() << "time spent in callback preparation: " << m.t_callback_prepare << " s." << endl;
    }
  }

  void Sqpmethod::printIteration(std::ostream &stream) const {
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
                                 double dx_norm, double rg, int ls_trials, bool ls_success) const {
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

  void Sqpmethod::reset_h(SqpmethodMemory& m) const {
    // Initial Hessian approximation of BFGS
    if (!exact_hessian_) {
      copy_n(B_init_.ptr(), Hsp_.nnz(), m.Bk);
    }
  }

  double Sqpmethod::getRegularization(const double* H) const {
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

  void Sqpmethod::regularize(double* H, double reg) const {
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


  void Sqpmethod::eval_h(SqpmethodMemory& m, const double* x, const double* lambda,
                         double sigma, double* H) const {
    try {
      calc_hess_l(m, x, m.p, &sigma, lambda, H);

      // Determing regularization parameter with Gershgorin theorem
      if (regularize_) {
        m.reg = getRegularization(H);
        if (m.reg > 0) {
          regularize(H, m.reg);
        }
      }

    } catch(exception& ex) {
      userOut<true, PL_WARN>() << "eval_h failed: " << ex.what() << endl;
      throw;
    }
  }

  void Sqpmethod::eval_g(SqpmethodMemory& m, const double* x, double* g) const {
    try {
      double time1 = clock();

      // Quick return if no constraints
      if (ng_==0) return;

      // Evaluate the function
      calc_g(m, x, m.p, g);

      double time2 = clock();
      m.t_eval_g += (time2-time1)/CLOCKS_PER_SEC;
      m.n_eval_g += 1;
    } catch(exception& ex) {
      userOut<true, PL_WARN>() << "eval_g failed: " << ex.what() << endl;
      throw;
    }
  }

  void Sqpmethod::eval_jac_g(SqpmethodMemory& m, const double* x, double* g, double* J) const {
    try {
      double time1 = clock();

      // Quich finish if no constraints
      if (ng_==0) return;

      // Evaluate the function
      calc_jac_g(m, x, m.p, g, J);

      double time2 = clock();
      m.t_eval_jac_g += (time2-time1)/CLOCKS_PER_SEC;
      m.n_eval_jac_g += 1;

    } catch(exception& ex) {
      userOut<true, PL_WARN>() << "eval_jac_g failed: " << ex.what() << endl;
      throw;
    }
  }

  void Sqpmethod::eval_grad_f(SqpmethodMemory& m, const double* x,
                              double* f, double* grad_f) const {
    try {
      double time1 = clock();

      // Evaluate the function
      calc_grad_f(m, x, m.p, f, grad_f);

      double time2 = clock();
      m.t_eval_grad_f += (time2-time1)/CLOCKS_PER_SEC;
      m.n_eval_grad_f += 1;

    } catch(exception& ex) {
      userOut<true, PL_WARN>() << "eval_grad_f failed: " << ex.what() << endl;
      throw;
    }
  }

  double Sqpmethod::eval_f(SqpmethodMemory& m, const double* x) const {
    try {
       // Log time
      double time1 = clock();

      // Evaluate the function
      double f;
      calc_f(m, x, m.p, &f);

      double time2 = clock();
      m.t_eval_f += (time2-time1)/CLOCKS_PER_SEC;
      m.n_eval_f += 1;

      return f;
    } catch(exception& ex) {
      userOut<true, PL_WARN>() << "eval_f failed: " << ex.what() << endl;
      throw;
    }
  }

  void Sqpmethod::solve_QP(SqpmethodMemory& m, const double* H, const double* g,
                           const double* lbx, const double* ubx,
                           const double* A, const double* lbA, const double* ubA,
                           double* x_opt, double* lambda_x_opt, double* lambda_A_opt) const {
    // Inputs
    fill_n(m.arg, qpsol_.n_in(), nullptr);
    m.arg[QPSOL_H] = H;
    m.arg[QPSOL_G] = g;
    m.arg[QPSOL_X0] = x_opt;
    m.arg[QPSOL_LBX] = lbx;
    m.arg[QPSOL_UBX] = ubx;
    m.arg[QPSOL_A] = A;
    m.arg[QPSOL_LBA] = lbA;
    m.arg[QPSOL_UBA] = ubA;

    // Outputs
    fill_n(m.res, qpsol_.n_out(), nullptr);
    m.res[QPSOL_X] = x_opt;
    m.res[QPSOL_LAM_X] = lambda_x_opt;
    m.res[QPSOL_LAM_A] = lambda_A_opt;

    // Solve the QP
    qpsol_(m.arg, m.res, m.iw, m.w, 0);
  }

  double Sqpmethod::
  primalInfeasibility(const double* x, const double* lbx, const double* ubx,
                      const double* g, const double* lbg, const double* ubg) const {
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
