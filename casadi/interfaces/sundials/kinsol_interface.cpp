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


#include "kinsol_interface.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_ROOTFINDER_KINSOL_EXPORT
  casadi_register_rootfinder_kinsol(Rootfinder::Plugin* plugin) {
    plugin->creator = KinsolInterface::creator;
    plugin->name = "kinsol";
    plugin->doc = KinsolInterface::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_ROOTFINDER_KINSOL_EXPORT casadi_load_rootfinder_kinsol() {
    Rootfinder::registerPlugin(casadi_register_rootfinder_kinsol);
  }

  /** \brief Kinsol solver class
   *
   * @copydoc Rootfinder_doc
   * You can provide an initial guess by setting output(0).\n
   * A good initial guess may be needed to avoid errors like
   * "The linear solver's setup function failed in an unrecoverable manner."
   *
   The constraints option expects an integer entry for each variable u:\n

   0 then no constraint is imposed on \p ui. \n
   1 then \p ui will be constrained to be \p ui >= 0.0. \n
   -1 then \p ui will be constrained to be \p ui <= 0.0. \n
   2 then \p ui will be constrained to be \p ui > 0.0. \n
   -2 then \p ui will be constrained to be \p ui < 0.0. \n

   *
   * \see Rootfinder for more information
   *
   */
  KinsolInterface::KinsolInterface(const std::string& name, const Function& f)
    : Rootfinder(name, f) {

    addOption("max_iter",                 OT_INTEGER, 0,
              "Maximum number of Newton iterations. Putting 0 sets the default value of KinSol.");
    addOption("abstol",                   OT_REAL, 1e-6, "Stopping criterion tolerance");
    addOption("linear_solver_type",       OT_STRING, "dense",
              "dense|banded|iterative|user_defined");
    addOption("upper_bandwidth",          OT_INTEGER);
    addOption("lower_bandwidth",          OT_INTEGER);
    addOption("max_krylov",               OT_INTEGER, 0);
    addOption("exact_jacobian",           OT_BOOLEAN, true);
    addOption("iterative_solver",         OT_STRING, "gmres", "gmres|bcgstab|tfqmr");
    addOption("f_scale",                  OT_REALVECTOR);
    addOption("u_scale",                  OT_REALVECTOR);
    addOption("pretype",                  OT_STRING, "none", "", "none|left|right|both");
    addOption("use_preconditioner",       OT_BOOLEAN, false); // precondition an iterative solver
    addOption("strategy",                 OT_STRING, "none", "Globalization strategy",
              "none|linesearch");
    addOption("disable_internal_warnings",   OT_BOOLEAN, false,
              "Disable KINSOL internal warning messages");
    addOption("monitor",      OT_STRINGVECTOR, GenericType(),  "", "eval_f|eval_djac", true);

    u_scale_ = 0;
    f_scale_ = 0;
    disable_internal_warnings_ = false;
  }

  KinsolInterface::~KinsolInterface() {
    if (u_scale_) N_VDestroy_Serial(u_scale_);
    if (f_scale_) N_VDestroy_Serial(f_scale_);
  }

  void KinsolInterface::init() {
    // Initialize the base classes
    Rootfinder::init();

    // Read options
    if (option("strategy")=="linesearch") {
      strategy_ = KIN_LINESEARCH;
    } else {
      casadi_assert(option("strategy")=="none");
      strategy_ = KIN_NONE;
    }

    // Use exact Jacobian?
    exact_jacobian_ = option("exact_jacobian");

    // Allocate N_Vectors
    if (u_scale_) N_VDestroy_Serial(u_scale_);
    if (f_scale_) N_VDestroy_Serial(f_scale_);
    u_scale_ = N_VNew_Serial(n_);
    f_scale_ = N_VNew_Serial(n_);

    // Set scaling factors on variables
    if (hasSetOption("u_scale")) {
      const vector<double>& u_scale = option("u_scale");
      casadi_assert(u_scale.size()==NV_LENGTH_S(u_scale_));
      copy(u_scale.begin(), u_scale.end(), NV_DATA_S(u_scale_));
    } else {
      N_VConst(1.0, u_scale_);
    }

    // Set scaling factors on equations
    if (hasSetOption("f_scale")) {
      const vector<double>& f_scale = option("f_scale");
      casadi_assert(f_scale.size()==NV_LENGTH_S(f_scale_));
      copy(f_scale.begin(), f_scale.end(), NV_DATA_S(f_scale_));
    } else {
      N_VConst(1.0, f_scale_);
    }

    // Disable internal warning messages?
    disable_internal_warnings_ = option("disable_internal_warnings");

    // Maximum number of iterations
    max_iter_ = option("max_iter");

    // Type of linear solver
    if (option("linear_solver_type")=="dense") {
      linear_solver_type_ = DENSE;
      if (exact_jacobian_) {
        // For storing Jacobian nonzeros
        alloc_w(jac_.nnz_out(0), true);
      }
    } else if (option("linear_solver_type")=="banded") {
      linear_solver_type_ = BANDED;
      upper_bandwidth_ = option("upper_bandwidth");
      lower_bandwidth_ = option("lower_bandwidth");
      if (exact_jacobian_) {
        // For storing Jacobian nonzeros
        alloc_w(jac_.nnz_out(0), true);
      }
    } else if (option("linear_solver_type")=="iterative") {
      linear_solver_type_ = ITERATIVE;
      maxl_ = option("max_krylov").toInt();
      if (option("iterative_solver")=="gmres") {
        iterative_solver_ = GMRES;
      } else if (option("iterative_solver")=="bcgstab") {
        iterative_solver_ = BCGSTAB;
      } else if (option("iterative_solver")=="tfqmr") {
        iterative_solver_ = TFQMR;
      } else {
        casadi_error("KINSOL: Unknown sparse solver");
      }
      if (exact_jacobian_) {
        // Form the Jacobian-times-vector function
        f_fwd_ = f_.derivative(1, 0);
        alloc(f_fwd_);
      }
      use_preconditioner_ = option("use_preconditioner");
      if (use_preconditioner_) {
        // Make sure that a Jacobian has been provided
        casadi_assert_message(!jac_.isNull(), "No Jacobian has been provided");

        // Make sure that a linear solver has been provided
        casadi_assert_message(!linsol_.isNull(), "No linear solver has been provided.");
      }
    } else if (option("linear_solver_type")=="user_defined") {
      linear_solver_type_ = USER_DEFINED;
      // Make sure that a Jacobian has been provided
      casadi_assert(!jac_.isNull());

      // Make sure that a linear solver has been provided
      casadi_assert(!linsol_.isNull());

      // Form the Jacobian-times-vector function
      f_fwd_ = f_.derivative(1, 0);
      alloc(f_fwd_);

      // Allocate space for Jacobian
      alloc_w(jac_.nnz_out(0), true);
    } else {
      casadi_error("Unknown linear solver");
    }

    // Stop criterion
    abstol_ = option("abstol");
  }

  void KinsolInterface::eval(const double** arg, double** res, int* iw, double* w, void* mem) {
    if (mem==0) {
      mem = alloc_mem();
      try {
        eval(arg, res, iw, w, mem);
      } catch (...) {
        free_mem(mem);
        throw;
      }
      free_mem(mem);
      return;
    }

    // Get memory block
    auto m = static_cast<KinsolMemory*>(mem);

    // Update IO references
    m->arg_ = arg;
    m->res_ = res;
    m->iw_ = iw;
    m->w_ = w;

    // Reset the counters
    m->t_func_ = 0;
    m->t_jac_ = 0;

    // Get the initial guess
    if (arg[iin_]) {
      copy(arg[iin_], arg[iin_]+nnz_in(iin_), NV_DATA_S(m->u_));
    } else {
      N_VConst(0.0, m->u_);
    }

    // Solve the nonlinear system of equations
    int flag = KINSol(m->mem_, m->u_, strategy_, u_scale_, f_scale_);
    if (flag<KIN_SUCCESS) kinsol_error("KINSol", flag);

    // Warn if not successful return
    if (verbose()) {
      if (flag!=KIN_SUCCESS) kinsol_error("KINSol", flag, false);
    }

    // Get the solution
    if (res[iout_]) {
      copy_n(NV_DATA_S(m->u_), nnz_out(iout_), res[iout_]);
    }

    // Evaluate auxiliary outputs
    if (n_out()>0) {
      // Temporary memory
      const double** arg1 = arg + n_in();
      double** res1 = res + n_out();

      // Evaluate f_
      copy_n(arg, n_in(), arg1);
      arg1[iin_] = NV_DATA_S(m->u_);
      copy_n(res, n_out(), res1);
      res1[iout_] = 0;
      f_(arg1, res1, iw, w, 0);
    }
  }

  void KinsolMemory::func(N_Vector u, N_Vector fval) {
    // Get time
    time1_ = clock();

    // Temporary memory
    const double** arg1 = arg_ + self.n_in();
    double** res1 = res_ + self.n_out();

    // Evaluate f_
    copy(arg_, arg_ + self.n_in(), arg1);
    arg1[self.iin_] = NV_DATA_S(u);
    fill_n(res1, self.n_out(), static_cast<double*>(0));
    res1[self.iout_] = NV_DATA_S(fval);
    self.f_(arg1, res1, iw_, w_, 0);

    // Print it, if requested
    if (self.monitored("eval_f")) {
      userOut() << "f = ";
      N_VPrint_Serial(fval);
    }

    // Make sure that all entries of the linear system are valid
    double *fdata = res1[self.iout_];
    for (int k=0; k<self.n_; ++k) {
      try {
        casadi_assert_message(!isnan(fdata[k]), "Nonzero " << k << " is not-a-number");
        casadi_assert_message(!isinf(fdata[k]), "Nonzero " << k << " is infinite");
      } catch(exception& ex) {
        stringstream ss;
        ss << ex.what() << endl;
        if (self.verbose()) {
          userOut() << "u = ";
          N_VPrint_Serial(u);

          // Print the expression for f[Jcol] if f is an SXFunction instance
          if (self.f_.is_a("sxfunction")) {
            vector<SX> res = self.f_(self.f_.sx_in());
            self.f_.print(ss);
            ss << "Equation " << k << " = " << res.at(0).at(k) << endl;
          }
        }

        throw CasadiException(ss.str());
      }
    }

    // Log time
    time2_ = clock();
    t_func_ += static_cast<double>(time2_-time1_)/CLOCKS_PER_SEC;
  }

  int KinsolMemory::func_wrapper(N_Vector u, N_Vector fval, void *user_data) {
    try {
      casadi_assert(user_data);
      auto this_ = static_cast<KinsolMemory*>(user_data);
      this_->func(u, fval);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "func failed: " << e.what() << endl;
      return 1;
    }
  }

  int KinsolMemory::djac_wrapper(long N, N_Vector u, N_Vector fu, DlsMat J,
                                            void *user_data, N_Vector tmp1, N_Vector tmp2) {
    try {
      casadi_assert(user_data);
      auto this_ = static_cast<KinsolMemory*>(user_data);
      this_->djac(N, u, fu, J, tmp1, tmp2);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "djac failed: " << e.what() << endl;;
      return 1;
    }
  }

  void KinsolMemory::djac(long N, N_Vector u, N_Vector fu, DlsMat J,
                                     N_Vector tmp1, N_Vector tmp2) {
    // Get time
    time1_ = clock();

    // Temporary memory
    const double** arg1 = arg_ + self.n_in();
    double** res1 = res_ + self.n_out();
    double* jac = w_ + self.jac_.sz_w();

    // Evaluate jac_
    copy(arg_, arg_ + self.n_in(), arg1);
    arg1[self.iin_] = NV_DATA_S(u);
    fill_n(res1, self.jac_.n_out(), static_cast<double*>(0));
    res1[0] = jac;
    self.jac_(arg1, res1, iw_, w_, 0);

    // Get sparsity and non-zero elements
    const int* colind = self.jac_.sparsity_out(0).colind();
    int ncol = self.jac_.size2_out(0);
    const int* row = self.jac_.sparsity_out(0).row();

    // Loop over columns
    for (int cc=0; cc<ncol; ++cc) {
      // Loop over non-zero entries
      for (int el=colind[cc]; el<colind[cc+1]; ++el) {
        // Get row
        int rr = row[el];

        // Set the element
        DENSE_ELEM(J, rr, cc) = jac[el];
      }
    }

    if (self.monitored("eval_djac")) {
      userOut() << "djac = ";
      PrintMat(J);
    }

    // Log time duration
    time2_ = clock();
    t_jac_ += static_cast<double>(time2_ - time1_)/CLOCKS_PER_SEC;
  }

  int KinsolMemory::bjac_wrapper(long N, long mupper, long mlower, N_Vector u,
                                            N_Vector fu, DlsMat J, void *user_data,
                                            N_Vector tmp1, N_Vector tmp2) {
    try {
      casadi_assert(user_data);
      auto this_ = static_cast<KinsolMemory*>(user_data);
      this_->bjac(N, mupper, mlower, u, fu, J, tmp1, tmp2);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "bjac failed: " << e.what() << endl;;
      return 1;
    }
  }

  void KinsolMemory::bjac(long N, long mupper, long mlower, N_Vector u, N_Vector fu,
                                     DlsMat J, N_Vector tmp1, N_Vector tmp2) {
    // Get time
    time1_ = clock();

    // Temporary memory
    const double** arg1 = arg_ + self.n_in();
    double** res1 = res_ + self.n_out();
    double* jac = w_ + self.jac_.sz_w();

    // Evaluate jac_
    copy(arg_, arg_ + self.jac_.n_in(), arg1);
    arg1[self.iin_] = NV_DATA_S(u);
    fill_n(res1, self.jac_.n_out(), static_cast<double*>(0));
    res1[0] = jac;
    self.jac_(arg1, res1, iw_, w_, 0);

    // Get sparsity and non-zero elements
    const int* colind = self.jac_.sparsity_out(0).colind();
    int ncol = self.jac_.size2_out(0);
    const int* row = self.jac_.sparsity_out(0).row();

    // Loop over cols
    for (int cc=0; cc<ncol; ++cc) {
      // Loop over non-zero entries
      for (int el=colind[cc]; el<colind[cc+1]; ++el) {
        // Get row
        int rr = row[el];

        // Set the element
        if (rr-cc>=-mupper && rr-cc<=mlower) {
          BAND_ELEM(J, rr, cc) = jac[el];
        }
      }
    }

    // Log time duration
    time2_ = clock();
    t_jac_ += static_cast<double>(time2_ - time1_)/CLOCKS_PER_SEC;
  }

  int KinsolMemory::jtimes_wrapper(N_Vector v, N_Vector Jv, N_Vector u, int* new_u,
                                     void *user_data) {
    try {
      casadi_assert(user_data);
      auto this_ = static_cast<KinsolMemory*>(user_data);
      this_->jtimes(v, Jv, u, new_u);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "jtimes failed: " << e.what() << endl;;
      return 1;
    }
  }

  void KinsolMemory::jtimes(N_Vector v, N_Vector Jv, N_Vector u, int* new_u) {
    // Get time
    time1_ = clock();

    // Temporary memory
    const double** arg1 = arg_ + self.n_in();
    double** res1 = res_ + self.n_out();

    // Evaluate f_fwd_
    copy(arg_, arg_ + self.n_in(), arg1);
    arg1[self.iin_] = NV_DATA_S(u);
    fill_n(arg1 + self.n_in(), self.n_in(), static_cast<const double*>(0));
    arg1[self.n_in()+self.iin_] = NV_DATA_S(v);
    fill_n(res1, self.f_fwd_.n_out(), static_cast<double*>(0));
    res1[self.n_out()] = NV_DATA_S(Jv);
    self.f_fwd_(arg1, res1, iw_, w_, 0);

    // Log time duration
    time2_ = clock();
    t_jac_ += static_cast<double>(time2_ - time1_)/CLOCKS_PER_SEC;
  }

  int KinsolMemory::
  psetup_wrapper(N_Vector u, N_Vector uscale, N_Vector fval, N_Vector fscale,
                 void* user_data, N_Vector tmp1, N_Vector tmp2) {
    try {
      casadi_assert(user_data);
      auto this_ = static_cast<KinsolMemory*>(user_data);
      this_->psetup(u, uscale, fval, fscale, tmp1, tmp2);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "psetup failed: " << e.what() << endl;;
      return 1;
    }
  }

  void KinsolMemory::
  psetup(N_Vector u, N_Vector uscale, N_Vector fval,
         N_Vector fscale, N_Vector tmp1, N_Vector tmp2) {
    // Get time
    time1_ = clock();

    // Temporary memory
    const double** arg1 = arg_ + self.n_in();
    double** res1 = res_ + self.n_out();
    double* jac = w_;
    double* w1 = w_ + self.jac_.nnz_out(0);

    // Evaluate jac_
    copy(arg_, arg_ + self.jac_.n_in(), arg1);
    arg1[self.iin_] = NV_DATA_S(u);
    fill_n(res1, self.jac_.n_out(), static_cast<double*>(0));
    res1[0] = jac;
    self.jac_(arg1, res1, iw_, w1, 0);

    // Get sparsity and non-zero elements
    const int* colind = self.jac_.sparsity_out(0).colind();
    int ncol = self.jac_.size2_out(0);
    const int* row = self.jac_.sparsity_out(0).row();

    // Log time duration
    time2_ = clock();
    t_lsetup_jac_ += static_cast<double>(time2_ - time1_)/CLOCKS_PER_SEC;

    // Prepare the solution of the linear system (e.g. factorize)
    fill_n(arg_, self.linsol_.n_in(), nullptr);
    fill_n(res_, self.linsol_.n_out(), nullptr);
    linsol_mem_ = Memory(self.linsol_, arg_, res_, iw_, w1, 0);
    self.linsol_.linsol_factorize(linsol_mem_, jac);

    // Log time duration
    time1_ = clock();
    t_lsetup_fac_ += static_cast<double>(time1_ - time2_)/CLOCKS_PER_SEC;
  }

  int KinsolMemory::psolve_wrapper(N_Vector u, N_Vector uscale, N_Vector fval,
                                     N_Vector fscale, N_Vector v, void* user_data, N_Vector tmp) {
    try {
      casadi_assert(user_data);
      auto this_ = static_cast<KinsolMemory*>(user_data);
      this_->psolve(u, uscale, fval, fscale, v, tmp);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "psolve failed: " << e.what() << endl;;
      return 1;
    }
  }

  void KinsolMemory::psolve(N_Vector u, N_Vector uscale, N_Vector fval,
                                       N_Vector fscale, N_Vector v, N_Vector tmp) {
    // Get time
    time1_ = clock();

    // Solve the factorized system
    self.linsol_.linsol_solve(linsol_mem_, NV_DATA_S(v));

    // Log time duration
    time2_ = clock();
    t_lsolve_ += static_cast<double>(time2_ - time1_)/CLOCKS_PER_SEC;
  }

  int KinsolMemory::lsetup_wrapper(KINMem kin_mem) {
    try {
      auto this_ = static_cast<KinsolMemory*>(kin_mem->kin_lmem);
      casadi_assert(this_);
      this_->lsetup(kin_mem);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "lsetup failed: " << e.what() << endl;;
      return -1;
    }
  }

  void KinsolMemory::lsetup(KINMem kin_mem) {
    N_Vector u =  kin_mem->kin_uu;
    N_Vector uscale = kin_mem->kin_uscale;
    N_Vector fval = kin_mem->kin_fval;
    N_Vector fscale = kin_mem->kin_fscale;
    N_Vector tmp1 = kin_mem->kin_vtemp1;
    N_Vector tmp2 = kin_mem->kin_vtemp2;
    psetup(u, uscale, fval, fscale, tmp1, tmp2);
  }

  int KinsolMemory::
  lsolve_wrapper(KINMem kin_mem, N_Vector x, N_Vector b, double *res_norm) {
    try {
      auto this_ = static_cast<KinsolMemory*>(kin_mem->kin_lmem);
      casadi_assert(this_);
      this_->lsolve(kin_mem, x, b, res_norm);
      return 0;
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "lsolve failed: " << e.what() << endl;;
      return -1;
    }
  }


  void KinsolMemory::lsolve(KINMem kin_mem, N_Vector x, N_Vector b, double *res_norm) {
    // Get vectors
    N_Vector u =  kin_mem->kin_uu;
    N_Vector uscale = kin_mem->kin_uscale;
    N_Vector fval = kin_mem->kin_fval;
    N_Vector fscale = kin_mem->kin_fscale;
    N_Vector tmp1 = kin_mem->kin_vtemp1;
    N_Vector tmp2 = kin_mem->kin_vtemp2;

    // Solve the linear system
    N_VScale(1.0, b, x);
    psolve(u, uscale, fval, fscale, x, tmp1);

    // Calculate residual
    jtimes(x, tmp2, u, 0);

    // Calculate the error in residual norm
    N_VLinearSum(1.0, b, -1.0, tmp2, tmp1);
    *res_norm = sqrt(N_VDotProd(tmp1, tmp1));
  }

  void KinsolMemory::
  ehfun_wrapper(int error_code, const char *module, const char *function,
                char *msg, void *eh_data) {
    try {
      casadi_assert(eh_data);
      auto this_ = static_cast<KinsolMemory*>(eh_data);
      this_->ehfun(error_code, module, function, msg);
    } catch(exception& e) {
      userOut<true, PL_WARN>() << "ehfun failed: " << e.what() << endl;
    }
  }

  void KinsolMemory::ehfun(int error_code, const char *module, const char *function,
                                      char *msg) {
    if (!self.disable_internal_warnings_) {
      userOut<true, PL_WARN>() << msg << endl;
    }
  }

  void KinsolInterface::kinsol_error(const string& module, int flag, bool fatal) {
    // Get the error message
    const char *id, *msg;
    switch (flag) {
    case KIN_SUCCESS:
      id = "KIN_SUCCES";
      msg = "KINSol succeeded; the scaled norm of F(u) is less than fnormtol";
    break;
  case KIN_INITIAL_GUESS_OK:
    id = "KIN_INITIAL_GUESS_OK";
    msg = "The guess u = u0 satisfied the system F(u) = 0 within the tolerances specified.";
    break;
    case KIN_STEP_LT_STPTOL:
      id = "KIN_STEP_LT_STPTOL";
      msg = "KINSol stopped based on scaled step length. This "
        "means that the current iterate may be an approximate solution of the "
        "given nonlinear system, but it is also quite possible that the algorithm"
        " is 'stalled' (making insufficient progress) near an invalid solution, "
        "or that the scalar scsteptol is too large.";
      break;
    case KIN_MEM_NULL:
      id = "KIN_MEM_NULL";
      msg = "The kinsol memory block pointer was NULL.";
      break;
    case KIN_ILL_INPUT:
      id = "KIN_ILL_INPUT";
      msg = "An input parameter was invalid.";
      break;
    case KIN_NO_MALLOC:
      id = "KIN_NO_MALLOC";
      msg = "The kinsol memory was not allocated by a call to KINCreate.";
      break;
    case KIN_LINESEARCH_NONCONV:
      id = "KIN_LINESEARCH_NONCONV";
      msg = "The line search algorithm was unable to find "
        "an iterate sufficiently distinct from the current iterate, or could not"
        " find an iterate satisfying the sufficient decrease condition. Failure"
        " to satisfy the sufficient decrease condition could mean the current "
        "iterate is 'close' to an approximate solution of the given nonlinear "
        "system, the difference approximation of the matrix-vector product J(u)v"
        " is inaccurate, or the real scalar scsteptol is too large.";
      break;
    case KIN_MAXITER_REACHED:
      id = "KIN_MAXITER_REACHED";
      msg = "The maximum number of nonlinear iterations "
        "has been reached.";
      break;
    case KIN_MXNEWT_5X_EXCEEDED:
      id = "KIN_MXNEWT_5X_EXCEEDED";
      msg = "Five consecutive steps have been taken that "
        "satisfy the inequality  || D_u p ||_L2 > 0.99 mxnewtstep, where p "
        "denotes the current step and mxnewtstep is a scalar upper bound on the "
        "scaled step length. Such a failure may mean that || D_F F(u)||_L2 "
        "asymptotes from above to a positive value, or the real scalar "
        "mxnewtstep is too small.";
      break;
    case KIN_LINESEARCH_BCFAIL:
      id = "KIN_LINESEARCH_BCFAIL";
      msg = "The line search algorithm was unable to satisfy "
        "the “beta-condition” for MXNBCF +1 nonlinear iterations (not necessarily "
        "consecutive), which may indicate the algorithm is making poor progress.";
      break;
    case KIN_LINSOLV_NO_RECOVERY:
      id = "KIN_LINSOLV_NO_RECOVERY";
      msg = "The user-supplied routine psolve encountered a"
        " recoverable error, but the preconditioner is already current.";
      break;
    case KIN_LINIT_FAIL:
      id = "KIN_LINIT_FAIL";
      msg = "The linear solver initialization routine (linit) encountered an error.";
    break;
    case KIN_LSETUP_FAIL:
      id = "KIN_LSETUP_FAIL";
      msg = "The user-supplied routine pset (used to set up the "
        "preconditioner data) encountered an unrecoverable error.";
      break;
    case KIN_LSOLVE_FAIL:
      id = "KIN_LSOLVE_FAIL";
      msg = "Either the user-supplied routine psolve "
        "(used to to solve the preconditioned linear system) encountered an "
        "unrecoverable error, or the linear solver routine (lsolve) "
        "encountered an error condition.";
      break;
    case KIN_SYSFUNC_FAIL:
      id = "KIN_SYSFUNC_FAIL";
      msg = "The system function failed in an unrecoverable manner.";
      break;
    case KIN_FIRST_SYSFUNC_ERR:
      id = "KIN_FIRST_SYSFUNC_ERR";
      msg = "The system function failed recoverably at the first call.";
      break;
    case KIN_REPTD_SYSFUNC_ERR:
      id = "KIN_REPTD_SYSFUNC_ERR";
      msg = "The system function had repeated recoverable errors. "
        "No recovery is possible.";
      break;
    default:
      id = "N/A";
      msg = 0;
    }

    // Construct message
    stringstream ss;
    if (msg==0) {
      ss << "Unknown " << (fatal? "error" : "warning") <<" (" << flag << ")"
        " from module \"" << module << "\".";
    } else {
      ss << "Module \"" << module << "\" returned flag \"" << id << "\"." << endl;
      ss << "The description of this flag is: " << endl;
      ss << "\"" << msg << "\"" << endl;
    }
    ss << "Consult KINSOL documentation for more information.";
    if (fatal) {
      casadi_error(ss.str())
        } else {
      casadi_warning(ss.str());
    }
  }

  KinsolMemory::KinsolMemory(KinsolInterface& s) : self(s) {
    // Current solution
    u_ = N_VNew_Serial(self.n_);

    // Create KINSOL memory block
    mem_ = KINCreate();

    // KINSOL bugfix
    KINMem kin_mem = KINMem(mem_);
    kin_mem->kin_inexact_ls = FALSE;

    // Set optional inputs
    int flag = KINSetUserData(mem_, this);
    casadi_assert_message(flag==KIN_SUCCESS, "KINSetUserData");

    // Set error handler function
    flag = KINSetErrHandlerFn(mem_, ehfun_wrapper, this);
    casadi_assert_message(flag==KIN_SUCCESS, "KINSetErrHandlerFn");

    // Initialize KINSOL
    flag = KINInit(mem_, func_wrapper, u_);
    casadi_assert(flag==KIN_SUCCESS);

    // Setting maximum number of Newton iterations
    flag = KINSetMaxNewtonStep(mem_, self.max_iter_);
    casadi_assert(flag==KIN_SUCCESS);

    // Set constraints
    if (!self.u_c_.empty()) {
      N_Vector domain  = N_VNew_Serial(self.n_);
      copy(self.u_c_.begin(), self.u_c_.end(), NV_DATA_S(domain));

      // Pass to KINSOL
      flag = KINSetConstraints(mem_, domain);
      casadi_assert(flag==KIN_SUCCESS);

      // Free the temporary vector
      N_VDestroy_Serial(domain);
    }

    switch (self.linear_solver_type_) {
    case KinsolInterface::DENSE:
      // Dense Jacobian
      flag = KINDense(mem_, self.n_);
      casadi_assert_message(flag==KIN_SUCCESS, "KINDense");

      if (self.exact_jacobian_) {
        flag = KINDlsSetDenseJacFn(mem_, djac_wrapper);
        casadi_assert_message(flag==KIN_SUCCESS, "KINDlsSetDenseJacFn");
      }
      break;
    case KinsolInterface::BANDED:
      // Banded Jacobian
      flag = KINBand(mem_, self.n_, self.upper_bandwidth_, self.lower_bandwidth_);
      casadi_assert_message(flag==KIN_SUCCESS, "KINBand");

      if (self.exact_jacobian_) {
        flag = KINDlsSetBandJacFn(mem_, bjac_wrapper);
        casadi_assert_message(flag==KIN_SUCCESS, "KINDlsBandJacFn");
      }
      break;
    case KinsolInterface::ITERATIVE:
      // Attach the sparse solver
      switch (self.iterative_solver_) {
      case KinsolInterface::GMRES:
        flag = KINSpgmr(mem_, self.maxl_);
        casadi_assert_message(flag==KIN_SUCCESS, "KINSpgmr");
        break;
      case KinsolInterface::BCGSTAB:
        flag = KINSpbcg(mem_, self.maxl_);
        casadi_assert_message(flag==KIN_SUCCESS, "KINSpbcg");
        break;
      case KinsolInterface::TFQMR:
        flag = KINSptfqmr(mem_, self.maxl_);
        casadi_assert_message(flag==KIN_SUCCESS, "KINSptfqmr");
        break;
      }

      // Attach functions for Jacobian information
      if (self.exact_jacobian_) {
        flag = KINSpilsSetJacTimesVecFn(mem_, jtimes_wrapper);
        casadi_assert_message(flag==KIN_SUCCESS, "KINSpilsSetJacTimesVecFn");
      }

      // Add a preconditioner
      if (self.use_preconditioner_) {
        flag = KINSpilsSetPreconditioner(mem_, psetup_wrapper, psolve_wrapper);
        casadi_assert(flag==KIN_SUCCESS);
      }
      break;
    case KinsolInterface::USER_DEFINED:
      // Set fields in the IDA memory
      KINMem kin_mem = KINMem(mem_);
      kin_mem->kin_lmem   = this;
      kin_mem->kin_lsetup = lsetup_wrapper;
      kin_mem->kin_lsolve = lsolve_wrapper;
      kin_mem->kin_setupNonNull = TRUE;
      break;
    }

    // Set stop criterion
    flag = KINSetFuncNormTol(mem_, self.abstol_);
    casadi_assert(flag==KIN_SUCCESS);
  }

  KinsolMemory::~KinsolMemory() {
    if (u_) N_VDestroy_Serial(u_);
    if (mem_) KINFree(&mem_);
  }

  void* KinsolInterface::alloc_mem() {
    return new KinsolMemory(*this);
  }

  void KinsolInterface::free_mem(void* mem) {
    delete static_cast<KinsolMemory*>(mem);
  }

} // namespace casadi
