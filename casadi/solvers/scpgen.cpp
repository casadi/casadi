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


#include "scpgen.hpp"
#include "casadi/core/core.hpp"
#include <ctime>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cfloat>
#include <cstdlib>
#ifdef HAVE_MKSTEMPS
#include <unistd.h>
#endif

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_NLPSOL_SCPGEN_EXPORT
      casadi_register_nlpsol_scpgen(Nlpsol::Plugin* plugin) {
    plugin->creator = Scpgen::creator;
    plugin->name = "scpgen";
    plugin->doc = Scpgen::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOL_SCPGEN_EXPORT casadi_load_nlpsol_scpgen() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_scpgen);
  }

  Scpgen::Scpgen(const std::string& name, const XProblem& nlp) : Nlpsol(name, nlp) {
    casadi_warning("SCPgen is under development");
    addOption("qpsol",         OT_STRING,   GenericType(),
              "The QP solver to be used by the SQP method");
    addOption("qpsol_options", OT_DICT, GenericType(),
              "Options to be passed to the QP solver");
    addOption("hessian_approximation", OT_STRING, "exact",
              "gauss-newton|exact");
    addOption("max_iter",          OT_INTEGER,      50,
              "Maximum number of SQP iterations");
    addOption("max_iter_ls",       OT_INTEGER,       1,
              "Maximum number of linesearch iterations");
    addOption("tol_pr",            OT_REAL,       1e-6,
              "Stopping criterion for primal infeasibility");
    addOption("tol_du",            OT_REAL,       1e-6,
              "Stopping criterion for dual infeasability");
    addOption("tol_reg",           OT_REAL,       1e-11,
              "Stopping criterion for regularization");
    addOption("tol_pr_step",       OT_REAL,       1e-6,
              "Stopping criterion for the step size");
    addOption("c1",                OT_REAL,       1e-4,
              "Armijo condition, coefficient of decrease in merit");
    addOption("beta",              OT_REAL,       0.8,
              "Line-search parameter, restoration factor of stepsize");
    addOption("merit_memsize",     OT_INTEGER,      4,
              "Size of memory to store history of merit function values");
    addOption("merit_start",       OT_REAL,      1e-8,
              "Lower bound for the merit function parameter");
    addOption("lbfgs_memory",      OT_INTEGER,     10,
              "Size of L-BFGS memory.");
    addOption("regularize",        OT_BOOLEAN,  false,
              "Automatic regularization of Lagrange Hessian.");
    addOption("print_header",      OT_BOOLEAN,   true,
              "Print the header with problem statistics");
    addOption("codegen",           OT_BOOLEAN,  false,
              "C-code generation");
    addOption("reg_threshold",     OT_REAL,      1e-8,
              "Threshold for the regularization.");
    addOption("name_x",      OT_STRINGVECTOR,  GenericType(),
              "Names of the variables.");
    addOption("print_x",           OT_INTEGERVECTOR,  GenericType(),
              "Which variables to print.");
    addOption("print_time",        OT_BOOLEAN, true,
              "Print information about execution time");

    // Monitors
    addOption("monitor",      OT_STRINGVECTOR, GenericType(),  "",
              "eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h|qp|dx", true);
  }


  Scpgen::~Scpgen() {
  }

  void Scpgen::init() {
    // Call the init method of the base class
    Nlpsol::init();

    // Read options
    max_iter_ = option("max_iter");
    max_iter_ls_ = option("max_iter_ls");
    c1_ = option("c1");
    beta_ = option("beta");
    lbfgs_memory_ = option("lbfgs_memory");
    tol_pr_ = option("tol_pr");
    tol_du_ = option("tol_du");
    tol_reg_ = option("tol_reg");
    regularize_ = option("regularize");
    codegen_ = option("codegen");
    reg_threshold_ = option("reg_threshold");
    print_time_ = option("print_time");
    tol_pr_step_ = option("tol_pr_step");
    merit_memsize_ = option("merit_memsize");
    merit_start_ = option("merit_start");
    string compiler = option("compiler");
    gauss_newton_ = option("hessian_approximation") == "gauss-newton";
    if (gauss_newton_) {
      casadi_assert(nlp_.output(NL_F).nnz()>1);
    } else {
      casadi_assert(nlp_.output(NL_F).nnz()==1);
    }

    // Name the components
    if (hasSetOption("name_x")) {
      name_x_ = option("name_x");
      casadi_assert(name_x_.size()==nx_);
    } else {
      stringstream ss;
      name_x_.resize(nx_);
      for (int i=0; i<nx_; ++i) {
        ss.str(string());
        ss << "x" << i;
        name_x_[i] = ss.str();
      }
    }

    // Components to print
    if (hasSetOption("print_x")) {
      print_x_ = option("print_x");
    } else {
      print_x_.resize(0);
    }

    Function fg = nlp_;
    if (!fg.is_a("mx_function")) {
      vector<MX> nlp_in = nlp_.mx_in();
      vector<MX> nlp_out = nlp_(nlp_in);
      fg = Function("fg", nlp_in, nlp_out);
    }

    // Generate lifting functions
    Function vdef_fcn, vinit_fcn;
    fg.generate_lifted(vdef_fcn, vinit_fcn);
    vinit_fcn_ = vinit_fcn;
    alloc(vinit_fcn_);

    // Extract the expressions
    vector<MX> vdef_in = vdef_fcn.mx_in();
    vector<MX> vdef_out = vdef_fcn(vdef_in);

    // Get the dimensions
    MX x = vdef_in.at(0);
    MX p = vdef_in.at(1);
    v_.resize(vdef_in.size()-2);
    for (int i=0; i<v_.size(); ++i) {
      v_[i].v = vdef_in.at(i+2);
      v_[i].v_def = vdef_out.at(i+2);
      v_[i].n = v_[i].v.nnz();
    }

    // Allocate memory, nonlfted problem
    alloc_w(nx_, true); // xk_
    alloc_w(ng_, true); // gk_
    alloc_w(nx_, true); // xk_
    alloc_w(ng_, true); // xg_
    alloc_w(nx_, true); // dxk_
    alloc_w(nx_, true); // lam_xk_
    alloc_w(nx_, true); // dlam_xk_
    alloc_w(ng_, true); // lam_gk_
    alloc_w(ng_, true); // dlam_gk_
    alloc_w(nx_, true); // gfk_
    alloc_w(nx_, true); // gL_

    // Allocate memory, lifted problem
    vm_.resize(v_.size());
    for (int i=0; i<v_.size(); ++i) {
      vm_[i].n = v_[i].n;
    }

    // Legacy
    qpH_times_du_.resize(nx_);

    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      it->init.resize(it->n, 0);
      it->opt.resize(it->n, 0);
      it->step.resize(it->n, 0);
      if (!gauss_newton_) {
        it->lam.resize(it->n, 0);
        it->dlam.resize(it->n, 0);
      }
    }

    // Line-search memory
    merit_mem_.resize(merit_memsize_);

    // Scalar objective function
    MX f;

    // Multipliers
    MX g_lam;

    // Definition of the lifted dual variables
    MX p_defL, gL_defL;

    if (gauss_newton_) {
      // Least square objective
      f = inner_prod(vdef_out[0], vdef_out[0])/2;
      gL_defL = vdef_out[0];
      ngn_ = gL_defL.nnz();
      alloc_w(ngn_, true); // b_gn_
      work_.resize(ngn_);
    } else {
      // Scalar objective function
      f = vdef_out[0];

      // Lagrange multipliers corresponding to the definition of the dependent variables
      stringstream ss;
      int i=0;
      for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
        ss.str(string());
        ss << "lam_x" << i++;
        it->v_lam = MX::sym(ss.str(), it->v.sparsity());
      }

      // Lagrange multipliers for the nonlinear constraints
      g_lam = MX::sym("g_lam", ng_);

      if (verbose_) {
        userOut() << "Allocated intermediate variables." << endl;
      }

      // Adjoint sweep to get the definitions of the lifted dual variables
      // (Equation 3.8 in Albersmeyer2010)
      vector<vector<MX> > fseed, fsens, aseed(1), asens(1);
      aseed[0].push_back(1.0);
      aseed[0].push_back(g_lam);
      for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
        aseed[0].push_back(it->v_lam);
      }
      vdef_fcn.derivative(vdef_in, vdef_out, fseed, fsens, aseed, asens, true);
      i=0;

      gL_defL = asens[0].at(i++);
      if (gL_defL.isNull()) gL_defL = MX::zeros(x.sparsity()); // Needed?

      p_defL = asens[0].at(i++);
      if (p_defL.isNull()) p_defL = MX::zeros(p.sparsity()); // Needed?

      for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
        it->v_defL = asens[0].at(i++);
        if (it->v_defL.isNull()) {
          it->v_defL = MX::zeros(it->v.sparsity());
        }
      }

      if (verbose_) {
        userOut() << "Generated the gradient of the Lagrangian." << endl;
      }
    }

    // Residual function

    // Inputs
    vector<MX> res_fcn_in;
    int n=0;
    res_fcn_in.push_back(x);             res_x_ = n++;
    res_fcn_in.push_back(p);             res_p_ = n++;
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      res_fcn_in.push_back(it->v);        it->res_var = n++;
    }
    if (!gauss_newton_) {
      res_fcn_in.push_back(g_lam);        res_g_lam_ = n++;
      for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
        res_fcn_in.push_back(it->v_lam);  it->res_lam = n++;
      }
    }

    // Outputs
    vector<MX> res_fcn_out;
    n=0;
    res_fcn_out.push_back(f);                              res_f_ = n++;
    res_fcn_out.push_back(gL_defL);                        res_gl_ = n++;
    res_fcn_out.push_back(vdef_out[1]);                    res_g_ = n++;
    res_fcn_out.push_back(p_defL);                         res_p_d_ = n++;
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      res_fcn_out.push_back(it->v_def - it->v);             it->res_d = n++;
    }

    if (!gauss_newton_) {
      for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
        res_fcn_out.push_back(it->v_defL - it->v_lam);     it->res_lam_d = n++;
      }
    }

    // Generate function
    Function res_fcn("res", res_fcn_in, res_fcn_out);
    if (verbose_) {
      userOut() << "Generated residual function ( " << res_fcn.countNodes() << " nodes)." << endl;
    }

    // Declare difference vector d and substitute out p and v
    stringstream ss;
    int i=0;
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      ss.str(string());
      ss << "d" << i++;
      it->d = MX::sym(ss.str(), it->v.sparsity());
      it->d_def = it->v_def - it->d;
    }

    // Declare difference vector lam_d and substitute out lam
    if (!gauss_newton_) {
      int i=0;
      for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
        ss.str(string());
        ss << "d_lam" << i++;
        it->d_lam = MX::sym(ss.str(), it->v.sparsity());
        it->d_defL = it->v_defL - it->d_lam;
      }
    }

    // Variables to be substituted and their definitions
    vector<MX> svar, sdef;
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      svar.push_back(it->v);
      sdef.push_back(it->d_def);
    }
    if (!gauss_newton_) {
      for (vector<Var>::reverse_iterator it=v_.rbegin(); it!=v_.rend(); ++it) {
        svar.push_back(it->v_lam);
        sdef.push_back(it->d_defL);
      }
    }

    vector<MX> ex(4);
    ex[0] = f;
    ex[1] = vdef_out[1];
    ex[2] = gL_defL;
    ex[3] = p_defL;

    substituteInPlace(svar, sdef, ex, false);
    i=0;
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      it->d_def = sdef[i++];
    }
    if (!gauss_newton_) {
      for (vector<Var>::reverse_iterator it=v_.rbegin(); it!=v_.rend(); ++it) {
        it->d_defL = sdef[i++];
      }
    }

    MX f_z = ex[0];
    MX g_z = ex[1];
    MX gL_z = ex[2];
    MX p_z = ex[3];

    // Modified function inputs
    vector<MX> mfcn_in;
    n=0;
    mfcn_in.push_back(p);                               mod_p_ = n++;
    mfcn_in.push_back(x);                               mod_x_ = n++;
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      mfcn_in.push_back(it->d);                          it->mod_var = n++;
    }

    // Modified function outputs
    n=0;
    vector<MX> mfcn_out;
    mfcn_out.push_back(g_z);                             mod_g_ = n++;

    // Add multipliers to function inputs
    if (!gauss_newton_) {
      n = mfcn_in.size();
      mfcn_in.push_back(g_lam);                          mod_g_lam_ = n++;
      for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
        mfcn_in.push_back(it->d_lam);                    it->mod_lam = n++;
      }
    }

    // Add gradient of the Lagrangian
    n = mfcn_out.size();
    mfcn_out.push_back(f_z);                             mod_f_ = n++;
    mfcn_out.push_back(gL_z);                            mod_gl_ = n++;

    // Lagrangian gradient function
    Function lgrad("lgrad", mfcn_in, mfcn_out);

    // Jacobian of the constraints
    MX jac = MX::jac(lgrad, mod_x_, mod_g_);
    log("Formed Jacobian of the constraints.");

    // Hessian of the Lagrangian
    MX hes = MX::jac(lgrad, mod_x_, mod_gl_, false, !gauss_newton_);
    if (gauss_newton_) {
      log("Formed square root of Gauss-Newton Hessian.");
    } else {
      log("Formed Hessian of the Lagrangian.");
    }

    // Matrices in the reduced QP
    n=0;
    vector<MX> mat_out;
    mat_out.push_back(jac);                             mat_jac_ = n++;
    mat_out.push_back(hes);                             mat_hes_ = n++;
    Function mat_fcn("mfcn", mfcn_in, mat_out);

    // Definition of intermediate variables
    n = mfcn_out.size();
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      mfcn_out.push_back(it->d_def);         it->mod_def = n++;
      if (!gauss_newton_) {
        mfcn_out.push_back(it->d_defL);      it->mod_defL = n++;
      }
    }

    // Modifier function
    Function mfcn("mfcn", mfcn_in, mfcn_out);

    // Directional derivative of Z
    vector<vector<MX> > mfcn_fwdSeed(1, mfcn_in), mfcn_fwdSens(1, mfcn_out);
    vector<vector<MX> > mfcn_adjSeed,            mfcn_adjSens;

    // Linearization in the d-direction (see Equation (2.12) in Alberspeyer2010)
    fill(mfcn_fwdSeed[0].begin(), mfcn_fwdSeed[0].end(), MX());
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      mfcn_fwdSeed[0][it->mod_var] = it->d;
      if (!gauss_newton_) {
        mfcn_fwdSeed[0][it->mod_lam] = it->d_lam;
      }
    }
    mfcn.derivative(mfcn_in, mfcn_out, mfcn_fwdSeed, mfcn_fwdSens,
                        mfcn_adjSeed, mfcn_adjSens, true);

    // Vector(s) b in Lifted Newton
    MX b_gf = densify(mfcn_fwdSens[0][mod_gl_]);
    MX b_g = densify(mfcn_fwdSens[0][mod_g_]);

    // Tangent function
    vector<MX> vec_fcn_out;
    n=0;
    vec_fcn_out.push_back(b_gf);                              vec_gf_ = n++;
    vec_fcn_out.push_back(b_g);                               vec_g_ = n++;
    casadi_assert(n==vec_fcn_out.size());

    Function vec_fcn("vec_fcn", mfcn_in, vec_fcn_out);
    if (verbose_) {
      userOut() << "Generated linearization function ( " << vec_fcn.countNodes()
           << " nodes)." << endl;
    }

    // Expression a + A*du in Lifted Newton (Section 2.1 in Alberspeyer2010)
    MX du = MX::sym("du", nx_);   // Step in u
    MX g_dlam;               // Step lambda_g
    if (!gauss_newton_) {
      g_dlam = MX::sym("g_dlam", g_lam.sparsity());
    }

    // Interpret the Jacobian-vector multiplication as a forward directional derivative
    fill(mfcn_fwdSeed[0].begin(), mfcn_fwdSeed[0].end(), MX());
    mfcn_fwdSeed[0][mod_x_] = du;
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      mfcn_fwdSeed[0][it->mod_var] = -it->d;
    }
    if (!gauss_newton_) {
      mfcn_fwdSeed[0][mod_g_lam_] = g_dlam;
      for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
        mfcn_fwdSeed[0][it->mod_lam] = -it->d_lam;
      }
    }

    mfcn.derivative(mfcn_in, mfcn_out, mfcn_fwdSeed, mfcn_fwdSens,
                        mfcn_adjSeed, mfcn_adjSens, true);

    // Step expansion function inputs
    n = mfcn_in.size();
    mfcn_in.push_back(du);                                 mod_du_ = n++;
    if (!gauss_newton_) {
      mfcn_in.push_back(g_dlam);                           mod_dlam_g_ = n++;
    }

    // Step expansion function outputs
    vector<MX> exp_fcn_out;
    n=0;
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      exp_fcn_out.push_back(mfcn_fwdSens[0][it->mod_def]); it->exp_def = n++;
    }

    if (!gauss_newton_) {
      for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
        exp_fcn_out.push_back(mfcn_fwdSens[0][it->mod_defL]); it->exp_defL = n++;
      }
    }

    // Step expansion function
    Function exp_fcn("exp_fcn", mfcn_in, exp_fcn_out);
    if (verbose_) {
      userOut() << "Generated step expansion function ( " << exp_fcn.countNodes() << " nodes)."
           << endl;
    }

    // Generate c code and load as DLL
    if (codegen_) {
      // Codegen the functions
      CodeGenerator gen;
      gen.add(res_fcn, "res_fcn");
      gen.add(mat_fcn, "mat_fcn");
      gen.add(vec_fcn, "vec_fcn");
      gen.add(exp_fcn, "exp_fcn");

      // Name of temporary file
      string cname;
#ifdef HAVE_MKSTEMPS
      // Preferred solution
      char cname_array[] = "tmp_casadi_scpgen_XXXXXX.c";
      if (mkstemps(cname_array, 3) == -1) {
        casadi_error("Failed to create a temporary file name");
      }
      cname = cname_array;
#else
      // Fallback, may result in deprecation warnings
      char* cname_array = tempnam(0, "scp.c");
      cname = cname_array;
      free(cname_array);
#endif

      // Generate code
      if (verbose_) {
        userOut() << "Generating \"" << cname << "\""  << endl;
      }
      string name = cname.substr(0, cname.find_first_of('.'));
      gen.generate(name);

      // Complile and run
      if (verbose_) {
        userOut() << "Starting compilation"  << endl;
      }
      time_t time1 = time(0);
      compiler_ = Compiler(cname, compilerplugin_, jit_options_);
      time_t time2 = time(0);
      double comp_time = difftime(time2, time1);
      if (verbose_) {
        userOut() << "Compilation completed after " << comp_time << " s."  << endl;
      }

      // Load the generated code
      res_fcn_ = Function::external("res_fcn", compiler_);
      mat_fcn_ = Function::external("mat_fcn", compiler_);
      vec_fcn_ = Function::external("vec_fcn", compiler_);
      exp_fcn_ = Function::external("exp_fcn", compiler_);
    } else {
      mat_fcn_ = mat_fcn;
      res_fcn_ = res_fcn;
      vec_fcn_ = vec_fcn;
      exp_fcn_ = exp_fcn;
    }
    alloc(mat_fcn_);
    alloc(res_fcn_);
    alloc(vec_fcn_);
    alloc(exp_fcn_);

    // Allocate QP data
    Sparsity sp_B_obj = mat_fcn_.output(mat_hes_).sparsity();
    qpH_ = DMatrix::zeros(sp_B_obj.T().zz_mtimes(sp_B_obj));
    qpA_ = mat_fcn_.output(mat_jac_);
    qpB_.resize(ng_);

    // QP solver options
    Dict qpsol_options;
    if (hasSetOption("qpsol_options")) {
      qpsol_options = option("qpsol_options");
    }

    // Allocate a QP solver
    qpsol_ = Function::qpsol("qpsol", option("qpsol"),
                                     {{"h", qpH_.sparsity()}, {"a", qpA_.sparsity()}},
                                     qpsol_options);
    alloc(qpsol_);
    if (verbose_) {
      userOut() << "Allocated QP solver." << endl;
    }

    // Residual
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      it->res.resize(it->d.nnz(), 0);
    }

    if (!gauss_newton_) {
      for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
        it->resL.resize(it->d_lam.nnz(), 0);
      }
    }

    if (verbose_) {
      userOut() << "NLP preparation completed" << endl;
    }

    // Header
    if (static_cast<bool>(option("print_header"))) {
      userOut() << "-------------------------------------------" << endl;
      userOut() << "This is casadi::SCPgen." << endl;
      if (gauss_newton_) {
        userOut() << "Using Gauss-Newton Hessian" << endl;
      } else {
        userOut() << "Using exact Hessian" << endl;
      }

      // Count the total number of variables
      int n_lifted = 0;
      for (vector<Var>::const_iterator i=v_.begin(); i!=v_.end(); ++i) {
        n_lifted += i->n;
      }

      userOut()
        << endl
        << "Number of reduced variables:               " << setw(9) << nx_ << endl
        << "Number of reduced constraints:             " << setw(9) << ng_ << endl
        << "Number of lifted variables/constraints:    " << setw(9) << n_lifted << endl
        << "Number of parameters:                      " << setw(9) << np_ << endl
        << "Total number of variables:                 " << setw(9) << (nx_+n_lifted) << endl
        << "Total number of constraints:               " << setw(9) << (ng_+n_lifted) << endl
        << endl;

      userOut()
        << "Iteration options:" << endl
        << "{ \"max_iter\":" << max_iter_ << ", "
        << "\"max_iter_ls\":" << max_iter_ls_ << ", "
        << "\"c1\":" << c1_ << ", "
        << "\"beta\":" << beta_ << ", "
        << "\"merit_memsize\":" << merit_memsize_ << ", "
        << "\"merit_start\":" << merit_start_ << ", "
        << "\"regularize\":" << regularize_ << ", "
        << endl << "  "
        << "\"tol_pr\":" << tol_pr_ << ", "
        << "\"tol_du\":" << tol_du_ << ", "
        << "\"tol_reg\":" << tol_reg_ << ", "
        << "\"reg_threshold\":" << reg_threshold_ << "}" << endl
        << endl;
    }
  }

  void Scpgen::reset(void* mem, const double**& arg, double**& res, int*& iw, double*& w) {
    // Reset the base classes
    Nlpsol::reset(mem, arg, res, iw, w);

    // Get work vectors, nonlifted problem
    xk_ = w; w += nx_;
    gk_ = w; w += ng_;
    dxk_ = w; w += nx_;
    lam_xk_ = w; w += nx_;
    dlam_xk_ = w; w += nx_;
    lam_gk_ = w; w += ng_;
    dlam_gk_ = w; w += ng_;
    gfk_ = w; w += nx_;
    gL_ = w; w += nx_;
    if (gauss_newton_) {
      b_gn_ = w; w += ngn_;
    }
  }

  void Scpgen::solve(void* mem) {
    for (int i=0; i<NLPSOL_NUM_IN; ++i) {
      const double *v;
      switch (i) {
      case NLPSOL_X0: v = x0_; break;
      case NLPSOL_P: v = p_; break;
      case NLPSOL_LBX: v = lbx_; break;
      case NLPSOL_UBX: v = ubx_; break;
      case NLPSOL_LBG: v = lbg_; break;
      case NLPSOL_UBG: v = ubg_; break;
      case NLPSOL_LAM_X0: v = lam_x0_; break;
      case NLPSOL_LAM_G0: v = lam_g0_; break;
      default: casadi_assert(0);
      }
      if (v) {
        setInputNZ(v, i);
      } else {
        setInput(0., i);
      }
    }

    // Check the provided inputs
    checkInputs(mem);
    if (v_.size()>0) {
      // Initialize lifted variables using the generated function
      fill_n(arg_, vinit_fcn_.n_in(), nullptr);
      arg_[0] = x0_;
      arg_[1] = p_;
      fill_n(res_, vinit_fcn_.n_out(), nullptr);
      for (int i=0; i<v_.size(); ++i) {
        res_[i] = getPtr(v_[i].init);
      }
      vinit_fcn_(arg_, res_, iw_, w_, 0);
    }
    if (verbose_) {
      userOut() << "Passed initial guess" << endl;
    }

    // Reset dual guess
    casadi_fill(lam_gk_, ng_, 0.);
    casadi_fill(dlam_gk_, ng_, 0.);
    casadi_fill(lam_xk_, nx_, 0.);
    casadi_fill(dlam_xk_, nx_, 0.);
    if (!gauss_newton_) {
      for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
        fill(it->lam.begin(), it->lam.end(), 0);
        fill(it->dlam.begin(), it->dlam.end(), 0);
      }
    }

    // Objective value
    fk_ = numeric_limits<double>::quiet_NaN();

    // Reset line-search
    fill(merit_mem_.begin(), merit_mem_.end(), 0.0);
    merit_ind_ = 0;

    // Current guess for the primal solution
    casadi_copy(x0_, nx_, xk_);
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      copy(it->init.begin(), it->init.end(), it->opt.begin());
    }

    // Get current time and reset timers
    double time1 = clock();
    t_eval_mat_ = t_eval_res_ = t_eval_vec_ = t_eval_exp_ = t_solve_qp_ = 0;

    // Initial evaluation of the residual function
    eval_res();

    // Number of SQP iterations
    int iter = 0;

    // Reset last step-size
    pr_step_ = 0;
    du_step_ = 0;

    // Reset line-search
    int ls_iter = 0;
    bool ls_success = true;

    // Reset regularization
    reg_ = 0;

    // Reset iteration message
    iteration_note_ = string();

    // MAIN OPTIMZATION LOOP
    while (true) {

      // Evaluate the vectors in the condensed QP
      eval_vec();

      // Evaluate the matrices in the condensed QP
      eval_mat();

      // 1-norm of the primal infeasibility
      double pr_inf = primalInfeasibility();

      // 1-norm of the dual infeasibility
      double du_inf = dualInfeasibility();

      // Print header occasionally
      if (iter % 10 == 0) printIteration(userOut());

      // Printing information about the actual iterate
      printIteration(userOut(), iter, fk_, pr_inf, du_inf, reg_, ls_iter, ls_success);

      // Checking convergence criteria
      bool converged = pr_inf <= tol_pr_ && pr_step_ <= tol_pr_step_ && reg_ <= tol_reg_;
      converged = converged && du_inf <= tol_du_;
      if (converged) {
        userOut()
          << endl
          << "casadi::SCPgen: Convergence achieved after " << iter << " iterations." << endl;
        break;
      }

      if (iter >= max_iter_) {
        userOut() << endl;
        userOut() << "casadi::SCPgen: Maximum number of iterations reached." << endl;
        break;
      }

      // Check if not-a-number
      if (fk_!=fk_ || pr_step_ != pr_step_ || pr_inf != pr_inf) {
        userOut() << "casadi::SCPgen: Aborted, nan detected" << endl;
        break;
      }

      // Start a new iteration
      iter++;

      // Regularize the QP
      if (regularize_) {
        regularize();
      }

      // Solve the condensed QP
      solve_qp();

      // Expand the step
      eval_exp();

      // Line-search to take the step
      line_search(ls_iter, ls_success);
    }

    double time2 = clock();
    t_mainloop_ = (time2-time1)/CLOCKS_PER_SEC;

    // Store optimal value
    userOut() << "optimal cost = " << fk_ << endl;

    // Save results to outputs
    output(NLPSOL_F).setNZ(&fk_);
    output(NLPSOL_X).setNZ(xk_);
    output(NLPSOL_LAM_G).setNZ(lam_gk_);
    output(NLPSOL_LAM_X).setNZ(lam_xk_);
    output(NLPSOL_G).setNZ(gk_);

    // Write timers
    if (print_time_) {
      userOut() << endl;
      userOut() << "time spent in eval_mat:    " << setw(9) << t_eval_mat_ << " s." << endl;
      userOut() << "time spent in eval_res:    " << setw(9) << t_eval_res_ << " s." << endl;
      userOut() << "time spent in eval_vec:    " << setw(9) << t_eval_vec_ << " s." << endl;
      userOut() << "time spent in eval_exp:    " << setw(9) << t_eval_exp_ << " s." << endl;
      userOut() << "time spent in solve_qp:    " << setw(9) << t_solve_qp_ << " s." << endl;
      userOut() << "time spent in main loop:   " << setw(9) << t_mainloop_ << " s." << endl;
    }

    // Save statistics
    stats_["iter_count"] = iter;

    userOut() << endl;

    for (int i=0; i<NLPSOL_NUM_OUT; ++i) {
      double **v;
      switch (i) {
      case NLPSOL_X: v = &x_; break;
      case NLPSOL_F: v = &f_; break;
      case NLPSOL_G: v = &g_; break;
      case NLPSOL_LAM_X: v = &lam_x_; break;
      case NLPSOL_LAM_G: v = &lam_g_; break;
      case NLPSOL_LAM_P: v = &lam_p_; break;
      default: casadi_assert(0);
      }
      if (*v) getOutputNZ(*v, i);
    }
  }

  double Scpgen::primalInfeasibility() {
    // L1-norm of the primal infeasibility
    double pr_inf = 0;

    // Simple bounds
    for (int i=0; i<nx_; ++i) pr_inf +=  std::max(xk_[i]-ubx(i), 0.);
    for (int i=0; i<nx_; ++i) pr_inf +=  std::max(lbx(i)-xk_[i], 0.);

    // Lifted variables
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      for (int i=0; i<it->n; ++i) pr_inf += ::fabs(it->res[i]);
    }

    // Nonlinear bounds
    for (int i=0; i<ng_; ++i) pr_inf += std::max(gk_[i]-ubg(i), 0.);
    for (int i=0; i<ng_; ++i) pr_inf += std::max(lbg(i)-gk_[i], 0.);

    return pr_inf;
  }

  double Scpgen::dualInfeasibility() {

    // L1-norm of the dual infeasibility
    double du_inf = 0;

    // Lifted variables
    for (int i=0; i<nx_; ++i) du_inf += ::fabs(gL_[i]);

    return du_inf;
  }

  void Scpgen::printIteration(std::ostream &stream) {
    stream << setw(4)  << "iter";
    stream << setw(14) << "objective";
    stream << setw(11) << "inf_pr";
    stream << setw(11) << "inf_du";
    stream << setw(11) << "pr_step";
    stream << setw(11) << "du_step";
    stream << setw(8) << "lg(rg)";
    stream << setw(3) << "ls";
    stream << ' ';

    // Print variables
    for (vector<int>::const_iterator i=print_x_.begin(); i!=print_x_.end(); ++i) {
      stream << setw(9) << name_x_.at(*i);
    }

    stream << endl;
    stream.unsetf(std::ios::floatfield);
  }

  void Scpgen::printIteration(std::ostream &stream, int iter, double obj,
                                      double pr_inf, double du_inf, double rg, int ls_trials,
                                      bool ls_success) {
    stream << setw(4) << iter;
    stream << scientific;
    stream << setw(14) << setprecision(6) << obj;
    stream << setw(11) << setprecision(2) << pr_inf;
    stream << setw(11);
    stream << setprecision(2) << du_inf;
    stream << setw(11) << setprecision(2) << pr_step_;
    stream << setw(11);
    stream << setprecision(2) << du_step_;
    stream << fixed;
    if (rg>0) {
      stream << setw(8) << setprecision(2) << log10(rg);
    } else {
      stream << setw(8) << "-";
    }
    stream << setw(3) << ls_trials;
    stream << (ls_success ? ' ' : 'F');

    // Print variables
    for (vector<int>::const_iterator i=print_x_.begin(); i!=print_x_.end(); ++i) {
      stream << setw(9) << setprecision(4) << xk_[*i];
    }

    // Print note
    if (!iteration_note_.empty()) {
      stream << "   " << iteration_note_;
      iteration_note_ = string();
    }

    stream.unsetf(std::ios::floatfield);
    stream << endl;
  }

  void Scpgen::eval_mat() {
    // Get current time
    double time1 = clock();

    // Pass parameters
    mat_fcn_.setInputNZ(p_, mod_p_);

    // Pass primal step/variables
    mat_fcn_.setInputNZ(xk_, mod_x_);
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      mat_fcn_.setInputNZ(it->res, it->mod_var);
    }

    // Pass dual steps/variables
    if (!gauss_newton_) {
      mat_fcn_.setInputNZ(lam_gk_, mod_g_lam_);
      for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
        mat_fcn_.setInputNZ(it->resL, it->mod_lam);
      }
    }

    // Evaluate condensed Hessian
    mat_fcn_.evaluate();

    // Get the Jacobian
    mat_fcn_.getOutput(qpA_, mat_jac_);

    if (gauss_newton_) {
      // Gauss-Newton Hessian
      const DMatrix& B_obj =  mat_fcn_.output(mat_hes_);
      fill(qpH_->begin(), qpH_->end(), 0);
      casadi_mul(B_obj.ptr(), B_obj.sparsity(), B_obj.ptr(), B_obj.sparsity(),
                 qpH_.ptr(), qpH_.sparsity(), getPtr(work_), true);

      // Gradient of the objective in Gauss-Newton
      casadi_fill(gfk_, nx_, 0.);
      casadi_mv(B_obj.ptr(), B_obj.sparsity(), b_gn_, gfk_, true);
    } else {
      // Exact Hessian
      mat_fcn_.getOutput(qpH_, mat_hes_);
    }

    // Calculate the gradient of the lagrangian
    const vector<double> &qpA_data = qpA_.data();
    const int* qpA_colind = qpA_.colind();
    int qpA_ncol = qpA_.size2();
    const int* qpA_row = qpA_.row();
    for (int i=0; i<nx_; ++i)  gL_[i] = gfk_[i] + lam_xk_[i];
    for (int cc=0; cc<qpA_ncol; ++cc) {
      for (int el=qpA_colind[cc]; el<qpA_colind[cc+1]; ++el) {
        int rr = qpA_row[el];
        gL_[cc] += qpA_data[el]*lam_gk_[rr];
      }
    }

    double time2 = clock();
    t_eval_mat_ += (time2-time1)/CLOCKS_PER_SEC;
  }

  void Scpgen::eval_res() {
    // Get current time
    double time1 = clock();

    // Pass parameters
    res_fcn_.setInputNZ(p_, res_p_);

    // Pass primal variables to the residual function for initial evaluation
    res_fcn_.setInputNZ(xk_, res_x_);
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      res_fcn_.setInputNZ(it->opt, it->res_var);
    }

    // Pass dual variables to the residual function for initial evaluation
    if (!gauss_newton_) {
      res_fcn_.setInput(0.0, res_g_lam_);
      for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
        res_fcn_.setInputNZ(it->lam, it->res_lam);
      }
    }

    // Evaluate residual function
    res_fcn_.evaluate();

    // Get objective
    fk_ = res_fcn_.output(res_f_).toScalar();

    // Get objective gradient
    if (gauss_newton_) {
      res_fcn_.getOutputNZ(b_gn_, res_gl_);
    } else {
      res_fcn_.getOutputNZ(gfk_, res_gl_);
    }

    // Get constraints
    res_fcn_.getOutputNZ(gk_, res_g_);

    // Get residuals
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      res_fcn_.getOutputNZ(it->res,  it->res_d);
      if (!gauss_newton_) {
        res_fcn_.getOutputNZ(it->resL, it->res_lam_d);
      }
    }

    // Parameter sensitivities
    res_fcn_.getOutput(output(NLPSOL_LAM_P), res_p_d_);

    double time2 = clock();
    t_eval_res_ += (time2-time1)/CLOCKS_PER_SEC;
  }

  void Scpgen::eval_vec() {
    // Get current time
    double time1 = clock();

    // Pass current parameter guess
    vec_fcn_.setInputNZ(p_, mod_p_);

    // Pass primal step/variables
    vec_fcn_.setInputNZ(xk_, mod_x_);
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      vec_fcn_.setInputNZ(it->res, it->mod_var);
    }

    // Pass dual steps/variables
    if (!gauss_newton_) {
      vec_fcn_.setInput(0.0, mod_g_lam_);
      for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
        vec_fcn_.setInputNZ(it->resL, it->mod_lam);
      }
    }

    // Evaluate to get QP
    vec_fcn_.evaluate();

    // Linear offset in the reduced QP
    transform(gk_, gk_ + ng_, vec_fcn_.output(vec_g_)->begin(),
              qpB_.begin(), std::minus<double>());

    // Gradient of the objective in the reduced QP
    if (gauss_newton_) {
      casadi_axpy(ngn_, -1., vec_fcn_.output(vec_gf_).ptr(), b_gn_);
    } else {
      casadi_axpy(nx_, -1., vec_fcn_.output(vec_gf_).ptr(), gfk_);
    }

    double time2 = clock();
    t_eval_vec_ += (time2-time1)/CLOCKS_PER_SEC;
  }

  void Scpgen::regularize() {
    casadi_assert(nx_==2);

    // Regularization
    reg_ = 0;

    // Check the smallest eigenvalue of the Hessian
    double a = qpH_.elem(0, 0);
    double b = qpH_.elem(0, 1);
    double c = qpH_.elem(1, 0);
    double d = qpH_.elem(1, 1);

    // Make sure no not a numbers
    casadi_assert(a==a && b==b && c==c &&  d==d);

    // Make sure symmetric
    if (b!=c) {
      casadi_assert_warning(fabs(b-c)<1e-10, "Hessian is not symmetric: " << b << " != " << c);
      qpH_.elem(1, 0) = c = b;
    }

    double eig_smallest = (a+d)/2 - std::sqrt(4*b*c + (a-d)*(a-d))/2;
    if (eig_smallest<reg_threshold_) {
      // Regularization
      reg_ = reg_threshold_-eig_smallest;
      qpH_(0, 0) += reg_;
      qpH_(1, 1) += reg_;
    }
  }

  void Scpgen::solve_qp() {
    // Get current time
    double time1 = clock();

    // Solve the QP
    qpsol_.setInputNZ(qpH_.ptr(), QPSOL_H);
    qpsol_.setInputNZ(gfk_, QPSOL_G);
    qpsol_.setInputNZ(qpA_.ptr(), QPSOL_A);
    double* qp_lbx = qpsol_.input(QPSOL_LBX).ptr();
    double* qp_ubx = qpsol_.input(QPSOL_UBX).ptr();
    double* qp_lba = qpsol_.input(QPSOL_LBA).ptr();
    double* qp_uba = qpsol_.input(QPSOL_UBA).ptr();

    // Get bounds
    casadi_copy(lbx_, nx_, qp_lbx);
    casadi_copy(ubx_, nx_, qp_ubx);
    casadi_copy(lbg_, ng_, qp_lba);
    casadi_copy(ubg_, ng_, qp_uba);
    casadi_axpy(nx_, -1., xk_, qp_lbx);
    casadi_axpy(nx_, -1., xk_, qp_ubx);
    casadi_axpy(ng_, -1., getPtr(qpB_), qp_lba);
    casadi_axpy(ng_, -1., getPtr(qpB_), qp_uba);
    qpsol_.evaluate();

    // Condensed primal step
    const DMatrix& du = qpsol_.output(QPSOL_X);
    copy(du->begin(), du->end(), dxk_);

    // Condensed dual step (simple bounds)
    const DMatrix& lam_x_new = qpsol_.output(QPSOL_LAM_X);
    copy(lam_x_new->begin(), lam_x_new->end(), dlam_xk_);
    casadi_axpy(nx_, -1., lam_xk_, dlam_xk_);

    // Condensed dual step (nonlinear bounds)
    const DMatrix& lam_g_new = qpsol_.output(QPSOL_LAM_A);
    copy(lam_g_new->begin(), lam_g_new->end(), dlam_gk_);
    casadi_axpy(ng_, -1., lam_gk_, dlam_gk_);

    double time2 = clock();
    t_solve_qp_ += (time2-time1)/CLOCKS_PER_SEC;
  }

  void Scpgen::line_search(int& ls_iter, bool& ls_success) {
    // Make sure that we have a decent direction
    if (!gauss_newton_) {
      // Get the curvature in the step direction
      double gain = casadi_qform(qpH_.ptr(), qpH_.sparsity(), dxk_);
      if (gain < 0) {
        iteration_note_ = "Hessian indefinite in the search direction";
      }
    }

    // Calculate penalty parameter of merit function
    sigma_ = merit_start_;
    sigma_ = std::max(sigma_, 1.01*norm_inf(qpsol_.output(QPSOL_LAM_X).data()));
    sigma_ = std::max(sigma_, 1.01*norm_inf(qpsol_.output(QPSOL_LAM_A).data()));

    // Calculate L1-merit function in the actual iterate
    double l1_infeas = primalInfeasibility();

    // Right-hand side of Armijo condition
    double F_sens = 0;
    for (int i=0; i<nx_; ++i) F_sens += dxk_[i] * gfk_[i];
    double L1dir = F_sens - sigma_ * l1_infeas;
    double L1merit = fk_ + sigma_ * l1_infeas;

    // Storing the actual merit function value in a list
    merit_mem_[merit_ind_] = L1merit;
    ++merit_ind_ %= merit_memsize_;

    // Stepsize
    double t = 1.0, t_prev = 0.0;

    // Merit function value in candidate
    double L1merit_cand = 0;

    // Reset line-search counter, success marker
    ls_iter = 0;
    ls_success = false;

    // Line-search
    //log("Starting line-search");

    // Line-search loop
    while (true) {

      // Take the primal step
      for (int i=0; i<nx_; ++i) xk_[i] += (t-t_prev) * dxk_[i];
      for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
        for (int i=0; i<it->n; ++i) it->opt[i] += (t-t_prev) * it->step[i];
      }

      // Take the dual step
      for (int i=0; i<ng_; ++i)  lam_gk_[i] += (t-t_prev) * dlam_gk_[i];
      for (int i=0; i<nx_; ++i)  lam_xk_[i] += (t-t_prev) * dlam_xk_[i];
      if (!gauss_newton_) {
        for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
          for (int i=0; i<it->n; ++i) it->lam[i] += (t-t_prev) * it->dlam[i];
        }
      }

      // Evaluate residual function to get objective and constraints
      // (and residuals for the next iteration)
      eval_res();
      ls_iter++;

      // Calculating merit-function in candidate
      l1_infeas = primalInfeasibility();
      L1merit_cand = fk_ + sigma_ * l1_infeas;

      // Calculating maximal merit function value so far
      double meritmax = *max_element(merit_mem_.begin(), merit_mem_.end());
      if (L1merit_cand <= meritmax + t * c1_ * L1dir) {

        // Accepting candidate
        ls_success = true;
        //log("Line-search completed, candidate accepted");
        break;
      }

      // Line-search not successful, but we accept it.
      if (ls_iter == max_iter_ls_) {
        //log("Line-search completed, maximum number of iterations");
        break;
      }

      // Backtracking
      t_prev = t;
      t = beta_ * t;
    }

    // Calculate primal step-size
    pr_step_ = 0;
    for (size_t i=0; i<nx_; ++i) pr_step_ += fabs(dxk_[i]);
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      for (vector<double>::const_iterator i=it->step.begin(); i!=it->step.end(); ++i)
        pr_step_ += fabs(*i);
    }
    pr_step_ *= t;

    // Calculate the dual step-size
    du_step_ = casadi_asum(ng_, dlam_gk_) + casadi_asum(nx_, dlam_xk_);
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      for (vector<double>::const_iterator i=it->dlam.begin(); i!=it->dlam.end(); ++i)
        du_step_ += fabs(*i);
    }
    du_step_ *= t;
  }

  void Scpgen::eval_exp() {
    // Get current time
    double time1 = clock();

    // Pass current parameter guess
    exp_fcn_.setInputNZ(p_, mod_p_);

    // Pass primal step/variables
    exp_fcn_.setInputNZ(dxk_, mod_du_);
    exp_fcn_.setInputNZ(xk_, mod_x_);
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      exp_fcn_.setInputNZ(it->res, it->mod_var);
    }

    // Pass dual step/variables
    if (!gauss_newton_) {
      exp_fcn_.setInputNZ(dlam_gk_, mod_dlam_g_);
      exp_fcn_.setInputNZ(lam_gk_, mod_g_lam_);
      for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
        exp_fcn_.setInputNZ(it->resL, it->mod_lam);
      }
    }

    // Perform the step expansion
    exp_fcn_.evaluate();

    // Expanded primal step
    for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
      const DMatrix& dv = exp_fcn_.output(it->exp_def);
      copy(dv->begin(), dv->end(), it->step.begin());
    }

    // Expanded dual step
    if (!gauss_newton_) {
      for (vector<Var>::iterator it=v_.begin(); it!=v_.end(); ++it) {
        const DMatrix& dlam_v = exp_fcn_.output(it->exp_defL);
        copy(dlam_v->begin(), dlam_v->end(), it->dlam.begin());
      }
    }

    double time2 = clock();
    t_eval_exp_ += (time2-time1)/CLOCKS_PER_SEC;
  }

} // namespace casadi
