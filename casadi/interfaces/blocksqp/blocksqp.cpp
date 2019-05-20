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


#include "blocksqp.hpp"
#include "casadi/core/casadi_misc.hpp"
#include "casadi/core/conic.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_NLPSOL_BLOCKSQP_EXPORT
  casadi_register_nlpsol_blocksqp(Nlpsol::Plugin* plugin) {
    plugin->creator = Blocksqp::creator;
    plugin->name = "blocksqp";
    plugin->doc = Blocksqp::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &Blocksqp::options_;
    plugin->deserialize = &Blocksqp::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOL_BLOCKSQP_EXPORT casadi_load_nlpsol_blocksqp() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_blocksqp);
  }

  Blocksqp::Blocksqp(const std::string& name, const Function& nlp)
    : Nlpsol(name, nlp) {
  }


  Blocksqp::~Blocksqp() {
    clear_mem();
  }

  const Options Blocksqp::options_
  = {{&Nlpsol::options_},
     {{"qpsol",
       {OT_STRING,
        "The QP solver to be used by the SQP method"}},
      {"qpsol_options",
       {OT_DICT,
        "Options to be passed to the QP solver"}},
      {"linsol",
       {OT_STRING,
        "The linear solver to be used by the QP method"}},
      {"print_header",
       {OT_BOOL,
        "Print solver header at startup"}},
      {"print_iteration",
       {OT_BOOL,
        "Print SQP iterations"}},
      {"eps",
       {OT_DOUBLE,
        "Values smaller than this are regarded as numerically zero"}},
      {"opttol",
       {OT_DOUBLE,
        "Optimality tolerance"}},
      {"nlinfeastol",
       {OT_DOUBLE,
        "Nonlinear feasibility tolerance"}},
      {"schur",
       {OT_BOOL,
        "Use qpOASES Schur compliment approach"}},
      {"globalization",
       {OT_BOOL,
        "Enable globalization"}},
      {"restore_feas",
       {OT_BOOL,
        "Use feasibility restoration phase"}},
      {"max_line_search",
       {OT_INT,
        "Maximum number of steps in line search"}},
      {"max_consec_reduced_steps",
       {OT_INT,
        "Maximum number of consecutive reduced steps"}},
      {"max_consec_skipped_updates",
       {OT_INT,
        "Maximum number of consecutive skipped updates"}},
      {"max_iter",
       {OT_INT,
        "Maximum number of SQP iterations"}},
      {"warmstart",
       {OT_BOOL,
        "Use warmstarting"}},
      {"qp_init",
       {OT_BOOL,
        "Use warmstarting"}},
      {"max_it_qp",
       {OT_INT,
        "Maximum number of QP iterations per SQP iteration"}},
      {"block_hess",
       {OT_INT,
        "Blockwise Hessian approximation?"}},
      {"hess_scaling",
       {OT_INT,
        "Scaling strategy for Hessian approximation"}},
      {"fallback_scaling",
       {OT_INT,
        "If indefinite update is used, the type of fallback strategy"}},
      {"max_time_qp",
       {OT_DOUBLE,
        "Maximum number of time in seconds per QP solve per SQP iteration"}},
      {"ini_hess_diag",
       {OT_DOUBLE,
        "Initial Hessian guess: diagonal matrix diag(iniHessDiag)"}},
      {"col_eps",
       {OT_DOUBLE,
        "Epsilon for COL scaling strategy"}},
      {"col_tau1",
       {OT_DOUBLE,
        "tau1 for COL scaling strategy"}},
      {"col_tau2",
       {OT_DOUBLE,
        "tau2 for COL scaling strategy"}},
      {"hess_damp",
       {OT_INT,
        "Activate Powell damping for BFGS"}},
      {"hess_damp_fac",
       {OT_DOUBLE,
        "Damping factor for BFGS Powell modification"}},
      {"hess_update",
       {OT_INT,
        "Type of Hessian approximation"}},
      {"fallback_update",
       {OT_INT,
        "If indefinite update is used, the type of fallback strategy"}},
      {"hess_lim_mem",
       {OT_INT,
        "Full or limited memory"}},
      {"hess_memsize",
       {OT_INT,
        "Memory size for L-BFGS updates"}},
      {"which_second_derv",
       {OT_INT,
        "For which block should second derivatives be provided by the user"}},
      {"skip_first_globalization",
       {OT_BOOL,
        "No globalization strategy in first iteration"}},
      {"conv_strategy",
       {OT_INT,
        "Convexification strategy"}},
      {"max_conv_qp",
       {OT_INT,
        "How many additional QPs may be solved for convexification per iteration?"}},
      {"max_soc_iter",
       {OT_INT,
        "Maximum number of SOC line search iterations"}},
      {"gamma_theta",
       {OT_DOUBLE,
        "Filter line search parameter, cf. IPOPT paper"}},
      {"gamma_f",
       {OT_DOUBLE,
        "Filter line search parameter, cf. IPOPT paper"}},
      {"kappa_soc",
       {OT_DOUBLE,
        "Filter line search parameter, cf. IPOPT paper"}},
      {"kappa_f",
       {OT_DOUBLE,
        "Filter line search parameter, cf. IPOPT paper"}},
      {"theta_max",
       {OT_DOUBLE,
        "Filter line search parameter, cf. IPOPT paper"}},
      {"theta_min",
       {OT_DOUBLE,
        "Filter line search parameter, cf. IPOPT paper"}},
      {"delta",
       {OT_DOUBLE,
        "Filter line search parameter, cf. IPOPT paper"}},
      {"s_theta",
       {OT_DOUBLE,
        "Filter line search parameter, cf. IPOPT paper"}},
      {"s_f",
       {OT_DOUBLE,
        "Filter line search parameter, cf. IPOPT paper"}},
      {"kappa_minus",
       {OT_DOUBLE,
        "Filter line search parameter, cf. IPOPT paper"}},
      {"kappa_plus",
       {OT_DOUBLE,
        "Filter line search parameter, cf. IPOPT paper"}},
      {"kappa_plus_max",
       {OT_DOUBLE,
        "Filter line search parameter, cf. IPOPT paper"}},
      {"delta_h0",
       {OT_DOUBLE,
        "Filter line search parameter, cf. IPOPT paper"}},
      {"eta",
       {OT_DOUBLE,
        "Filter line search parameter, cf. IPOPT paper"}},
      {"obj_lo",
       {OT_DOUBLE,
        "Lower bound on objective function [-inf]"}},
      {"obj_up",
       {OT_DOUBLE,
        "Upper bound on objective function [inf]"}},
      {"rho",
       {OT_DOUBLE,
        "Feasibility restoration phase parameter"}},
      {"zeta",
       {OT_DOUBLE,
        "Feasibility restoration phase parameter"}},
      {"print_maxit_reached",
       {OT_BOOL,
        "Print error when maximum number of SQP iterations reached"}}
     }
  };

  void Blocksqp::init(const Dict& opts) {
    // Call the init method of the base class
    Nlpsol::init(opts);

    // Set default options
    //string qpsol_plugin = "qpoases";
    //Dict qpsol_options;
    linsol_plugin_ = "ma27";
    print_header_ = true;
    print_iteration_ = true;
    eps_ = 1.0e-16;
    opttol_ = 1.0e-6;
    nlinfeastol_ = 1.0e-6;
    schur_ = true;
    globalization_ = true;
    restore_feas_ = true;
    max_line_search_ = 20;
    max_consec_reduced_steps_ = 100;
    max_consec_skipped_updates_ = 100;
    max_iter_ = 100;
    warmstart_ = false;
    max_it_qp_ = 5000;
    block_hess_ = true;
    hess_scaling_ = 2;
    fallback_scaling_ = 4;
    max_time_qp_ = 10000.0;
    ini_hess_diag_ = 1.0;
    col_eps_ = 0.1;
    col_tau1_ = 0.5;
    col_tau2_ = 1.0e4;
    hess_damp_ = 1;
    hess_damp_fac_ = 0.2;
    hess_update_ = 1;
    fallback_update_ = 2;
    hess_lim_mem_ = 1;
    hess_memsize_ = 20;
    which_second_derv_ = 0;
    skip_first_globalization_ = false;
    conv_strategy_ = 0;
    max_conv_qp_ = 1;
    max_soc_iter_ = 3;
    gamma_theta_ = 1.0e-5;
    gamma_f_ = 1.0e-5;
    kappa_soc_ = 0.99;
    kappa_f_ = 0.999;
    theta_max_ = 1.0e7;
    theta_min_ = 1.0e-5;
    delta_ = 1.0;
    s_theta_ = 1.1;
    s_f_ = 2.3;
    kappa_minus_ = 0.333;
    kappa_plus_ = 8.0;
    kappa_plus_max_ = 100.0;
    delta_h0_ = 1.0e-4;
    eta_ = 1.0e-4;
    obj_lo_ = -inf;
    obj_up_ = inf;
    rho_ = 1.0e3;
    zeta_ = 1.0e-3;
    print_maxit_reached_ = true;
    qp_init_ = true;

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="qpsol") {
        //qpsol_plugin = op.second.to_string();
        casadi_warning("Option 'qpsol' currently not supported, ignored");
      } else if (op.first=="qpsol_options") {
        //qpsol_options = op.second;
        casadi_warning("Option 'qpsol_options' currently not supported, ignored");
      } else if (op.first=="linsol") {
        linsol_plugin_ = string(op.second);
      } else if (op.first=="print_header") {
        print_header_ = op.second;
      } else if (op.first=="print_iteration") {
          print_iteration_ = op.second;
      } else if (op.first=="eps") {
        eps_ = op.second;
      } else if (op.first=="opttol") {
        opttol_ = op.second;
      } else if (op.first=="nlinfeastol") {
        nlinfeastol_ = op.second;
      } else if (op.first=="schur") {
        schur_ = op.second;
      } else if (op.first=="globalization") {
        globalization_ = op.second;
      } else if (op.first=="restore_feas") {
        restore_feas_ = op.second;
      } else if (op.first=="max_line_search") {
        max_line_search_ = op.second;
      } else if (op.first=="max_consec_reduced_steps") {
        max_consec_reduced_steps_ = op.second;
      } else if (op.first=="max_consec_skipped_updates") {
        max_consec_skipped_updates_ = op.second;
      } else if (op.first=="max_iter") {
        max_iter_ = op.second;
      } else if (op.first=="warmstart") {
        warmstart_ = op.second;
      } else if (op.first=="qp_init") {
        qp_init_ = op.second;
      } else if (op.first=="max_it_qp") {
        max_it_qp_ = op.second;
      } else if (op.first=="block_hess") {
        block_hess_ = op.second;
      } else if (op.first=="hess_scaling") {
        hess_scaling_ = op.second;
      } else if (op.first=="fallback_scaling") {
        fallback_scaling_ = op.second;
      } else if (op.first=="max_time_qp") {
        max_time_qp_ = op.second;
      } else if (op.first=="ini_hess_diag") {
        ini_hess_diag_ = op.second;
      } else if (op.first=="col_eps") {
        col_eps_ = op.second;
      } else if (op.first=="col_tau1") {
        col_tau1_ = op.second;
      } else if (op.first=="col_tau2") {
        col_tau2_ = op.second;
      } else if (op.first=="hess_damp") {
        hess_damp_ = op.second;
      } else if (op.first=="hess_damp_fac") {
        hess_damp_fac_ = op.second;
      } else if (op.first=="hess_update") {
        hess_update_ = op.second;
      } else if (op.first=="fallback_update") {
        fallback_update_ = op.second;
      } else if (op.first=="hess_lim_mem") {
        hess_lim_mem_ = op.second;
      } else if (op.first=="hess_memsize") {
        hess_memsize_ = op.second;
      } else if (op.first=="which_second_derv") {
        which_second_derv_ = op.second;
      } else if (op.first=="skip_first_globalization") {
        skip_first_globalization_ = op.second;
      } else if (op.first=="conv_strategy") {
        conv_strategy_ = op.second;
      } else if (op.first=="max_conv_qp") {
        max_conv_qp_ = op.second;
      } else if (op.first=="max_soc_iter") {
        max_soc_iter_ = op.second;
      } else if (op.first=="gamma_theta") {
        gamma_theta_ = op.second;
      } else if (op.first=="gamma_f") {
        gamma_f_ = op.second;
      } else if (op.first=="kappa_soc") {
        kappa_soc_ = op.second;
      } else if (op.first=="kappa_f") {
        kappa_f_ = op.second;
      } else if (op.first=="theta_max") {
        theta_max_ = op.second;
      } else if (op.first=="theta_min") {
        theta_min_ = op.second;
      } else if (op.first=="delta") {
        delta_ = op.second;
      } else if (op.first=="s_theta") {
        s_theta_ = op.second;
      } else if (op.first=="s_f") {
        s_f_ = op.second;
      } else if (op.first=="kappa_minus") {
        kappa_minus_ = op.second;
      } else if (op.first=="kappa_plus") {
        kappa_plus_ = op.second;
      } else if (op.first=="kappa_plus_max") {
        kappa_plus_max_ = op.second;
      } else if (op.first=="delta_h0") {
        delta_h0_ = op.second;
      } else if (op.first=="eta") {
        eta_ = op.second;
      } else if (op.first=="obj_lo") {
        obj_lo_ = op.second;
      } else if (op.first=="obj_up") {
        obj_up_ = op.second;
      } else if (op.first=="rho") {
        rho_ = op.second;
      } else if (op.first=="zeta") {
        zeta_ = op.second;
      } else if (op.first=="print_maxit_reached") {
        print_maxit_reached_ = op.second;
      }
    }

    // If we compute second constraints derivatives switch to
    // finite differences Hessian (convenience)
    if (which_second_derv_ == 2) {
      hess_update_ = 4;
      block_hess_ = true;
    }

    // If we don't use limited memory BFGS we need to store only one vector.
    if (!hess_lim_mem_) hess_memsize_ = 1;
    if (!schur_ && hess_update_ == 1) {
      print("***WARNING: SR1 update only works with qpOASES Schur complement version. "
             "Using BFGS updates instead.\n");
      hess_update_ = 2;
      hess_scaling_ = fallback_scaling_;
    }

    // Setup feasibility restoration phase
    if (restore_feas_) {
      // get orignal nlp
      Function nlp = oracle_;
      vector<MX> resv;
      vector<MX> argv = nlp.mx_in();

      // build fesibility restoration phase nlp
      MX p = MX::sym("p", nlp.size1_in("x"));
      MX s = MX::sym("s", nlp.size1_out("g"));

      MX d = fmin(1.0, 1.0/abs(p)) * (argv.at(0) - p);
      MX f_rp = 0.5 * rho_ * dot(s, s) + zeta_/2.0 * dot(d, d);
      MX g_rp = nlp(argv).at(1) - s;

      MXDict nlp_rp = {{"x", MX::vertcat({argv.at(0), s})},
                       {"p", MX::vertcat({argv.at(1), p})},
                       {"f", f_rp},
                       {"g", g_rp}};

      // Set options for the SQP method for the restoration problem
      Dict solver_options;
      solver_options["globalization"] = true;
      solver_options["which_second_derv"] = 0;
      solver_options["restore_feas"] = false;
      solver_options["hess_update"] = 2;
      solver_options["hess_lim_mem"] = 1;
      solver_options["hess_scaling"] = 2;
      solver_options["opttol"] = opttol_;
      solver_options["nlinfeastol"] = nlinfeastol_;
      solver_options["max_iter"] = 1;
      solver_options["print_time"] = false;
      solver_options["print_header"] = false;
      solver_options["print_iteration"] = false;
      solver_options["print_maxit_reached"] = false;

      // Create and initialize solver for the restoration problem
      rp_solver_ = nlpsol("rpsolver", "blocksqp", nlp_rp, solver_options);
    }

    // Setup NLP functions
    create_function("nlp_fg", {"x", "p"}, {"f", "g"});
    Function gf_jg = create_function("nlp_gf_jg", {"x", "p"},
                                     {"f", "g", "grad:f:x", "jac:g:x"});
    Asp_ = gf_jg.sparsity_out("jac_g_x");

    if (!block_hess_) {
      // No block-structured Hessian
      blocks_ = {0, nx_};
      which_second_derv_ = 0;
      Hsp_ = Sparsity::dense(nx_, nx_);
    } else {
      // Detect block structure

      // Get the sparsity pattern for the Hessian of the Lagrangian
      Function grad_lag = oracle_.factory("grad_lag",
                                          {"x", "p", "lam:f", "lam:g"}, {"grad:gamma:x"},
                                          {{"gamma", {"f", "g"}}});
      Hsp_ = grad_lag.sparsity_jac("x", "grad_gamma_x", false, true);

      // Make sure diagonal exists
      Hsp_ = Hsp_ + Sparsity::diag(nx_);

      // Find the strongly connected components of the Hessian
      // Unlike Sparsity::scc, assume ordered
      const casadi_int* colind = Hsp_.colind();
      const casadi_int* row = Hsp_.row();
      blocks_.push_back(0);
      casadi_int ind = 0;
      while (ind < nx_) {
        // Find the next cutoff
        casadi_int next=ind+1;
        while (ind<next && ind<nx_) {
          for (casadi_int k=colind[ind]; k<colind[ind+1]; ++k) next = max(next, 1+row[k]);
          ind++;
        }
        blocks_.push_back(next);
      }
    }

    // Number of blocks
    nblocks_ = blocks_.size()-1;

    // Blocksizes
    dim_.resize(nblocks_);
    casadi_int max_size = 0;
    nnz_H_ = 0;
    for (casadi_int i=0; i<nblocks_; ++i) {
      dim_[i] = blocks_[i+1]-blocks_[i];
      max_size = max(max_size, dim_[i]);
      nnz_H_ += dim_[i]*dim_[i];
    }

    create_function("nlp_hess_l", {"x", "p", "lam:f", "lam:g"},
                    {"hess:gamma:x:x"}, {{"gamma", {"f", "g"}}});
    exact_hess_lag_sp_ = get_function("nlp_hess_l").sparsity_out(0);

    if (verbose_) casadi_message(str(nblocks_) + " blocks of max size " + str(max_size) + ".");

    // Allocate a QP solver
    //casadi_assert(!qpsol_plugin.empty(), "'qpsol' option has not been set");
    //qpsol_ = conic("qpsol", qpsol_plugin, {{"h", Hsp_}, {"a", Asp_}},
    //               qpsol_options);
    //alloc(qpsol_);

    // Allocate memory
    alloc_w(Asp_.nnz(), true); // jac
    alloc_w(nx_, true); // lam_xk
    alloc_w(ng_, true); // lam_gk
    alloc_w(ng_, true); // gk
    alloc_w(nx_, true); // grad_fk
    alloc_w(nx_, true); // grad_lagk
    alloc_w(nx_+ng_, true); // lam_qp
    alloc_w(nblocks_, true); // delta_norm
    alloc_w(nblocks_, true); // delta_norm_old
    alloc_w(nblocks_, true); // delta_gamma
    alloc_w(nblocks_, true); // delta_gamma_old
    alloc_w(nblocks_, true); // delta_h
    alloc_w(nx_, true); // trial_xk
    alloc_w(nx_, true); // lbx_qp
    alloc_w(nx_, true); // ubx_qp
    alloc_w(ng_, true); // lba_qp
    alloc_w(ng_, true); // uba_qp
    alloc_w(ng_, true); // jac_times_dxk
    alloc_w(nx_*hess_memsize_, true); // deltaMat
    alloc_w(nx_*hess_memsize_, true); // gammaMat
    alloc_w(Asp_.nnz(), true); // jac_g
    alloc_w(nnz_H_, true); // hess_lag
    alloc_iw(nblocks_, true); // noUpdateCounter

    // Allocate block diagonal Hessian(s)
    casadi_int n_hess = hess_update_==1 || hess_update_==4 ? 2 : 1;
    alloc_res(nblocks_*n_hess, true);
    alloc_w(n_hess*nnz_H_, true);
    alloc_iw(nnz_H_ + (nx_+1) + nx_, true); // hessIndRow
    alloc_w(exact_hess_lag_sp_.nnz(), true); // exact_hess_lag
  }

  int Blocksqp::init_mem(void* mem) const {
    if (Nlpsol::init_mem(mem)) return 1;
    auto m = static_cast<BlocksqpMemory*>(mem);

    // Create qpOASES memory
    if (schur_) {
      m->qpoases_mem = new QpoasesMemory();
      m->qpoases_mem->linsol_plugin = linsol_plugin_;
    }

    m->qp = nullptr;
    m->colind.resize(Asp_.size2()+1);
    m->row.resize(Asp_.nnz());
    return 0;
  }

  void Blocksqp::set_work(void* mem, const double**& arg, double**& res,
                                   casadi_int*& iw, double*& w) const {
    auto m = static_cast<BlocksqpMemory*>(mem);

    // Set work in base classes
    Nlpsol::set_work(mem, arg, res, iw, w);

    // Temporary memory
    m->jac = w; w += Asp_.nnz();
    m->lam_xk = w; w += nx_;
    m->lam_gk = w; w += ng_;
    m->gk = w; w += ng_;
    m->grad_fk = w; w += nx_;
    m->grad_lagk = w; w += nx_;
    m->lam_qp = w; w += nx_+ng_;
    m->delta_norm = w; w += nblocks_;
    m->delta_norm_old = w; w += nblocks_;
    m->delta_gamma = w; w += nblocks_;
    m->delta_gamma_old = w; w += nblocks_;
    m->delta_h = w; w += nblocks_;
    m->trial_xk = w; w += nx_;
    m->lbx_qp = w; w += nx_;
    m->ubx_qp = w; w += nx_;
    m->lba_qp = w; w += ng_;
    m->uba_qp = w; w += ng_;
    m->jac_times_dxk = w; w += ng_;
    m->deltaMat = w; w += nx_*hess_memsize_;
    m->gammaMat = w; w += nx_*hess_memsize_;
    m->jac_g = w; w += Asp_.nnz();
    m->hess_lag = w; w += nnz_H_;
    m->hessIndRow = reinterpret_cast<int*>(iw); iw += nnz_H_ + (nx_+1) + nx_;
    m->noUpdateCounter = iw; iw += nblocks_;

    // First Hessian
    m->hess1 = res; res += nblocks_;
    for (casadi_int b=0; b<nblocks_; b++) {
      m->hess1[b] = w; w += dim_[b]*dim_[b];
    }

    // Second Hessian, for SR1 or finite differences
    if (hess_update_ == 1 || hess_update_ == 4) {
      m->hess2 = res; res += nblocks_;
      for (casadi_int b=0; b<nblocks_; b++) {
        m->hess2[b] = w; w += dim_[b]*dim_[b];
      }
    } else {
      m->hess2 = nullptr;
    }

    m->exact_hess_lag = w; w += exact_hess_lag_sp_.nnz();
  }

  int Blocksqp::solve(void* mem) const {
    auto m = static_cast<BlocksqpMemory*>(mem);
    auto d_nlp = &m->d_nlp;

    casadi_int ret = 0;

    // Create problem evaluation object
    vector<casadi_int> blocks = blocks_;

    /*-------------------------------------------------*/
    /* Create blockSQP method object and run algorithm */
    /*-------------------------------------------------*/
    m->itCount = 0;
    m->qpItTotal = 0;
    m->qpIterations = 0;
    m->qpIterations2 = 0;
    m->qpResolve = 0;
    m->rejectedSR1 = 0;
    m->hessSkipped = 0;
    m->hessDamped = 0;
    m->averageSizingFactor = 0.0;
    m->nFunCalls = 0;
    m->nDerCalls = 0;
    m->nRestHeurCalls = 0;
    m->nRestPhaseCalls = 0;

    m->nTotalUpdates = 0;
    m->nTotalSkippedUpdates = 0;

    casadi_int maxblocksize = 1;

    for (casadi_int k=0; k<nblocks_+1; k++) {
      if (k > 0)
        if (blocks_[k] - blocks_[k-1] > maxblocksize)
          maxblocksize = blocks_[k] - blocks_[k-1];
    }

    if (hess_lim_mem_ && hess_memsize_ == 0)
      const_cast<Blocksqp*>(this)->hess_memsize_ = maxblocksize;

    // Reset the SQP metod
    reset_sqp(m);

    // Free existing memory, if any
    if (qp_init_) {
      delete m->qp;
      m->qp = nullptr;
    }
    if (!m->qp) {
      if (schur_) {
        m->qp = new qpOASES::SQProblemSchur(nx_, ng_, qpOASES::HST_UNKNOWN, 50,
                                            m->qpoases_mem,
                                            QpoasesInterface::qpoases_init,
                                            QpoasesInterface::qpoases_sfact,
                                            QpoasesInterface::qpoases_nfact,
                                            QpoasesInterface::qpoases_solve);
      } else {
        m->qp = new qpOASES::SQProblem(nx_, ng_);
      }
    }

    // Print header and information about the algorithmic parameters
    if (print_header_) printInfo(m);

    // Open output files
    initStats(m);
    initIterate(m);

    // Initialize filter with pair (maxConstrViolation, objLowerBound)
    initializeFilter(m);

    // Primal-dual initial guess
    casadi_copy(d_nlp->lam, nx_, m->lam_xk);
    casadi_scal(nx_, -1., m->lam_xk);
    casadi_copy(d_nlp->lam + nx_, ng_, m->lam_gk);
    casadi_scal(ng_, -1., m->lam_gk);

    casadi_copy(m->lam_xk, nx_, m->lam_qp);
    casadi_copy(m->lam_gk, ng_, m->lam_qp+nx_);

    ret = run(m, max_iter_, warmstart_);

    m->success = ret==0;
    m->ret_ = ret;

    if (ret==1) {
      if (print_maxit_reached_) print("***WARNING: Maximum number of iterations reached\n");
      m->unified_return_status = SOLVER_RET_LIMITED;
    }


    // Get optimal cost
    d_nlp->f = m->obj;
    // Get constraints at solution
    casadi_copy(m->gk, ng_, d_nlp->z + nx_);
    // Get dual solution (simple bounds)
    if (d_nlp->lam) {
      casadi_copy(m->lam_xk, nx_, d_nlp->lam);
      casadi_scal(nx_, -1., d_nlp->lam);
    }
    // Get dual solution (nonlinear bounds)
    casadi_copy(m->lam_gk, ng_, d_nlp->lam + nx_);
    casadi_scal(ng_, -1., d_nlp->lam + nx_);
    return 0;
  }

  casadi_int Blocksqp::run(BlocksqpMemory* m, casadi_int maxIt, casadi_int warmStart) const {
    casadi_int it, infoQP = 0;
    bool skipLineSearch = false;
    bool hasConverged = false;

    if (warmStart == 0 || m->itCount == 0) {
      // SQP iteration 0

      /// Set initial Hessian approximation
      calcInitialHessian(m);

      /// Evaluate all functions and gradients for xk_0
      switch (evaluate(m, &m->obj, m->gk, m->grad_fk, m->jac_g)) {
        case -1:
          m->unified_return_status = SOLVER_RET_NAN;
          return -1;
        case 0:
          break;
        default:
          return 1;
      }
      m->nDerCalls++;

      /// Check if converged
      hasConverged = calcOptTol(m);
      if (print_iteration_) printProgress(m);
      updateStats(m);
      if (hasConverged) {
        if (print_iteration_ && m->steptype < 2) {
          print("\n***CONVERGENCE ACHIEVED!***\n");
        }
        return 0;
      }
      m->itCount++;
    }

    /*
     * SQP Loop: during first iteration, m->itCount = 1
     */
    for (it=0; it<maxIt; it++) {
      /// Solve QP subproblem with qpOASES or QPOPT
      updateStepBounds(m, 0);
      infoQP = solveQP(m, m->dxk, m->lam_qp);

      if (infoQP == 1) {
          // 1.) Maximum number of iterations reached
          print("***WARNING: Maximum number of QP iterations exceeded.***\n");
      } else if (infoQP == 2 || infoQP > 3) {
          // 2.) QP error (e.g., unbounded), solve again with pos.def. diagonal matrix (identity)
          print("***WARNING: QP error. Solve again with identity matrix.***\n");
          resetHessian(m);
          infoQP = solveQP(m, m->dxk, m->lam_qp);
          if (infoQP) {
            // If there is still an error, terminate.
            print("***WARNING: QP error. Stop.***\n");
            return -1;
          } else {
            m->steptype = 1;
          }
        } else if (infoQP == 3) {
        // 3.) QP infeasible, try to restore feasibility
        bool qpError = true;
        skipLineSearch = true; // don't do line search with restoration step

        // Try to reduce constraint violation by heuristic
        if (m->steptype < 2) {
          print("***WARNING: QP infeasible. Trying to reduce constraint violation ...");
          qpError = feasibilityRestorationHeuristic(m);
          if (!qpError) {
            m->steptype = 2;
            print("Success.***\n");
          } else {
            print("Failure.***\n");
          }
        }

        // Invoke feasibility restoration phase
        //if (qpError && m->steptype < 3 && restore_feas_)
        if (qpError && restore_feas_ && m->cNorm > 0.01 * nlinfeastol_) {
          print("***Start feasibility restoration phase.***\n");
          m->steptype = 3;
          qpError = feasibilityRestorationPhase(m);
        }

        // If everything failed, abort.
        if (qpError) {
          print("***WARNING: QP error. Stop.\n");
          return -1;
        }
      }

      /// Determine steplength alpha
      if (!globalization_ || (skip_first_globalization_ && m->itCount == 1)) {
        // No globalization strategy, but reduce step if function cannot be evaluated
        if (fullstep(m)) {
          print("***WARNING: Constraint or objective could "
                         "not be evaluated at new point. Stop.***\n");
          return -1;
        }
        m->steptype = 0;
      } else if (globalization_ && !skipLineSearch) {
        // Filter line search based on Waechter et al., 2006 (Ipopt paper)
        if (filterLineSearch(m) || m->reducedStepCount > max_consec_reduced_steps_) {
          // Filter line search did not produce a step. Now there are a few things we can try ...
          bool lsError = true;

          // Heuristic 1: Check if the full step reduces the KKT error by at
          // least kappaF, if so, accept the step.
          lsError = kktErrorReduction(m);
          if (!lsError)
            m->steptype = -1;

          // Heuristic 2: Try to reduce constraint violation by closing
          // continuity gaps to produce an admissable iterate
          if (lsError && m->cNorm > 0.01 * nlinfeastol_ && m->steptype < 2) {
            // Don't do this twice in a row!
            print("***WARNING: Steplength too short. "
                  "Trying to reduce constraint violation...");

            // Integration over whole time interval
            lsError = feasibilityRestorationHeuristic(m);
            if (!lsError) {
                m->steptype = 2;
                print("Success.***\n");
              } else {
              print("***WARNING: Failed.***\n");
            }
          }

          // Heuristic 3: Recompute step with a diagonal Hessian
          if (lsError && m->steptype != 1 && m->steptype != 2) {
            // After closing continuity gaps, we already take a step with initial Hessian.
            // If this step is not accepted then this will cause an infinite loop!

            print("***WARNING: Steplength too short. "
                  "Trying to find a new step with identity Hessian.***\n");
            m->steptype = 1;

            resetHessian(m);
            continue;
          }

          // If this does not yield a successful step, start restoration phase
          if (lsError && m->cNorm > 0.01 * nlinfeastol_ && restore_feas_) {
            print("***WARNING: Steplength too short. "
                           "Start feasibility restoration phase.***\n");
            m->steptype = 3;

            // Solve NLP with minimum norm objective
            lsError = feasibilityRestorationPhase(m);
          }

          // If everything failed, abort.
          if (lsError) {
            print("***WARNING: Line search error. Stop.***\n");
            return -1;
          }
        } else {
          m->steptype = 0;
        }
      }

      /// Calculate "old" Lagrange gradient: gamma = dL(xi_k, lambda_k+1)
      calcLagrangeGradient(m, m->gamma, 0);

      /// Evaluate functions and gradients at the new xk
      (void)evaluate(m, &m->obj, m->gk, m->grad_fk, m->jac_g);
      m->nDerCalls++;

      /// Check if converged
      hasConverged = calcOptTol(m);

      /// Print one line of output for the current iteration
      if (print_iteration_) printProgress(m);
      updateStats(m);
      if (hasConverged && m->steptype < 2) {
        if (print_iteration_) print("\n***CONVERGENCE ACHIEVED!***\n");
        m->itCount++;
        return 0; //Convergence achieved!
      }

      /// Calculate difference of old and new Lagrange gradient:
      // gamma = -gamma + dL(xi_k+1, lambda_k+1)
      calcLagrangeGradient(m, m->gamma, 1);

      /// Revise Hessian approximation
      if (hess_update_ < 4 && !hess_lim_mem_) {
        calcHessianUpdate(m, hess_update_, hess_scaling_);
      } else if (hess_update_ < 4 && hess_lim_mem_) {
        calcHessianUpdateLimitedMemory(m, hess_update_, hess_scaling_);
      } else if (hess_update_ == 4) {
        calcHessianUpdateExact(m);
      }

      // If limited memory updates  are used, set pointers deltaXi and
      // gamma to the next column in deltaMat and gammaMat
      updateDeltaGamma(m);

      m->itCount++;
      skipLineSearch = false;
    }

    return 1;
  }

  /**
   * Compute gradient of Lagrangian or difference of Lagrangian gradients (sparse version)
   *
   * flag == 0: output dL(xk, lambda)
   * flag == 1: output dL(xi_k+1, lambda_k+1) - L(xi_k, lambda_k+1)
   * flag == 2: output dL(xi_k+1, lambda_k+1) - df(xi)
   */
  void Blocksqp::
  calcLagrangeGradient(BlocksqpMemory* m,
    const double* lam_x, const double* lam_g,
    const double* grad_f, const double *jacNz,
    double* grad_lag, casadi_int flag) const {

    // Objective gradient
    if (flag == 0) {
      casadi_copy(grad_f, nx_, grad_lag);
    } else if (flag == 1) {
      casadi_scal(nx_, -1., grad_lag);
      casadi_axpy(nx_, 1., grad_f, grad_lag);
    } else {
      casadi_fill(grad_lag, nx_, 0.);
    }

    // - lambdaT * constrJac
    const casadi_int* jacIndRow = Asp_.row();
    const casadi_int* jacIndCol = Asp_.colind();
    for (casadi_int iVar=0; iVar<nx_; iVar++) {
      for (casadi_int iCon=jacIndCol[iVar]; iCon<jacIndCol[iVar+1]; iCon++) {
        grad_lag[iVar] -= lam_g[jacIndRow[iCon]] * jacNz[iCon];
      }
    }

    // - lambdaT * simpleBounds
    casadi_axpy(nx_, -1., lam_x, grad_lag);
  }

  /**
   * Wrapper if called with standard arguments
   */
  void Blocksqp::
  calcLagrangeGradient(BlocksqpMemory* m, double* grad_lag, casadi_int flag) const {
    calcLagrangeGradient(m, m->lam_xk, m->lam_gk, m->grad_fk, m->jac_g,
      grad_lag, flag);
  }


  /**
   * Compute optimality conditions:
   * ||grad_lag(xk,lambda)||_infty / (1 + ||lambda||_infty) <= TOL
   * and
   * ||constrViolation||_infty / (1 + ||xi||_infty) <= TOL
   */
  bool Blocksqp::calcOptTol(BlocksqpMemory* m) const {
    auto d_nlp = &m->d_nlp;
    // scaled norm of Lagrangian gradient
    calcLagrangeGradient(m, m->grad_lagk, 0);
    m->gradNorm = casadi_norm_inf(nx_, m->grad_lagk);
    m->tol = m->gradNorm/(1.0+fmax(casadi_norm_inf(nx_, m->lam_xk),
                                   casadi_norm_inf(ng_, m->lam_gk)));

    // norm of constraint violation
    m->cNorm  = lInfConstraintNorm(m, d_nlp->z, m->gk);
    m->cNormS = m->cNorm /(1.0 + casadi_norm_inf(nx_, d_nlp->z));

    return m->tol <= opttol_ && m->cNormS <= nlinfeastol_;
  }

  void Blocksqp::printInfo(BlocksqpMemory* m) const {
    char hessString1[100];
    char hessString2[100];
    char globString[100];
    char qpString[100];

    /* QP Solver */
    if (schur_)
      strcpy(qpString, "sparse, Schur complement approach");
    else
      strcpy(qpString, "sparse, reduced Hessian factorization");

    /* Globalization */
    if (globalization_) {
      strcpy(globString, "filter line search");
    } else {
      strcpy(globString, "none (full step)");
    }

    /* Hessian approximation */
    if (block_hess_ && (hess_update_ == 1 || hess_update_ == 2))
      strcpy(hessString1, "block ");
    else
      strcpy(hessString1, "");

    if (hess_lim_mem_ && (hess_update_ == 1 || hess_update_ == 2))
      strcat(hessString1, "L-");

    /* Fallback Hessian */
    if (hess_update_ == 1 || hess_update_ == 4
      || (hess_update_ == 2 && !hess_damp_)) {
        strcpy(hessString2, hessString1);

        /* Fallback Hessian update type */
        if (fallback_update_ == 0) {
          strcat(hessString2, "Id");
        } else if (fallback_update_ == 1) {
          strcat(hessString2, "SR1");
        } else if (fallback_update_ == 2) {
          strcat(hessString2, "BFGS");
        } else if (fallback_update_ == 4) {
          strcat(hessString2, "Finite differences");
        }

        /* Fallback Hessian scaling */
        if (fallback_scaling_ == 1) {
          strcat(hessString2, ", SP");
        } else if (fallback_scaling_ == 2) {
          strcat(hessString2, ", OL");
        } else if (fallback_scaling_ == 3) {
          strcat(hessString2, ", mean");
        } else if (fallback_scaling_ == 4) {
          strcat(hessString2, ", selective sizing");
        }
      } else {
        strcpy(hessString2, "-");
    }

    /* First Hessian update type */
    if (hess_update_ == 0) {
      strcat(hessString1, "Id");
    } else if (hess_update_ == 1) {
      strcat(hessString1, "SR1");
    } else if (hess_update_ == 2) {
      strcat(hessString1, "BFGS");
    } else if (hess_update_ == 4) {
      strcat(hessString1, "Finite differences");
    }

    /* First Hessian scaling */
    if (hess_scaling_ == 1) {
      strcat(hessString1, ", SP");
    } else if (hess_scaling_ == 2) {
      strcat(hessString1, ", OL");
    } else if (hess_scaling_ == 3) {
      strcat(hessString1, ", mean");
    } else if (hess_scaling_ == 4) {
      strcat(hessString1, ", selective sizing");
    }

    print("\n+---------------------------------------------------------------+\n");
    print("| Starting blockSQP with the following algorithmic settings:    |\n");
    print("+---------------------------------------------------------------+\n");
    print("| qpOASES flavor            | %-34s|\n", qpString);
    print("| Globalization             | %-34s|\n", globString);
    print("| 1st Hessian approximation | %-34s|\n", hessString1);
    print("| 2nd Hessian approximation | %-34s|\n", hessString2);
    print("+---------------------------------------------------------------+\n\n");
  }

  void Blocksqp::
  acceptStep(BlocksqpMemory* m, double alpha) const {
    acceptStep(m, m->dxk, m->lam_qp, alpha, 0);
  }

  void Blocksqp::
  acceptStep(BlocksqpMemory* m, const double* deltaXi,
    const double* lambdaQP, double alpha, casadi_int nSOCS) const {
    auto d_nlp = &m->d_nlp;
    double lStpNorm;

    // Current alpha
    m->alpha = alpha;
    m->nSOCS = nSOCS;

    // Set new x by accepting the current trial step
    for (casadi_int k=0; k<nx_; k++) {
      d_nlp->z[k] = m->trial_xk[k];
      m->dxk[k] = alpha * deltaXi[k];
    }

    // Store the infinity norm of the multiplier step
    m->lambdaStepNorm = 0.0;
    for (casadi_int k=0; k<nx_; k++)
      if ((lStpNorm = fabs(alpha*lambdaQP[k] - alpha*m->lam_xk[k])) > m->lambdaStepNorm)
        m->lambdaStepNorm = lStpNorm;
    for (casadi_int k=0; k<ng_; k++)
      if ((lStpNorm = fabs(alpha*lambdaQP[nx_+k] - alpha*m->lam_gk[k])) > m->lambdaStepNorm)
        m->lambdaStepNorm = lStpNorm;

    // Set new multipliers
    for (casadi_int k=0; k<nx_; k++)
      m->lam_xk[k] = (1.0 - alpha)*m->lam_xk[k] + alpha*lambdaQP[k];
    for (casadi_int k=0; k<ng_; k++)
      m->lam_gk[k] = (1.0 - alpha)*m->lam_gk[k] + alpha*lambdaQP[nx_+k];

    // Count consecutive reduced steps
    if (m->alpha < 1.0)
      m->reducedStepCount++;
    else
      m->reducedStepCount = 0;
  }

  void Blocksqp::
  reduceStepsize(BlocksqpMemory* m, double *alpha) const {
    *alpha = (*alpha) * 0.5;
  }

  void Blocksqp::
  reduceSOCStepsize(BlocksqpMemory* m, double *alphaSOC) const {
    auto d_nlp = &m->d_nlp;
    // Update bounds on linearized constraints for the next SOC QP:
    // That is different from the update for the first SOC QP!
    for (casadi_int i=0; i<ng_; i++) {
      double lbg = d_nlp->lbz[i + nx_];
      double ubg = d_nlp->ubz[i + nx_];
      if (lbg != inf) {
        m->lba_qp[i] = *alphaSOC * m->lba_qp[i] - m->gk[i];
      } else {
        m->lba_qp[i] = inf;
      }

      if (ubg != inf) {
        m->uba_qp[i] = *alphaSOC * m->uba_qp[i] - m->gk[i];
      } else {
        m->uba_qp[i] = inf;
      }
    }

    *alphaSOC *= 0.5;
  }

  /**
   * Take a full Quasi-Newton step, except when integrator fails:
   * xk = xk + deltaXi
   * lambda = lambdaQP
   */
  casadi_int Blocksqp::fullstep(BlocksqpMemory* m) const {
    auto d_nlp = &m->d_nlp;
    double alpha;
    double objTrial, cNormTrial;

    // Backtracking line search
    alpha = 1.0;
    for (casadi_int k=0; k<10; k++) {
      // Compute new trial point
      for (casadi_int i=0; i<nx_; i++)
        m->trial_xk[i] = d_nlp->z[i] + alpha * m->dxk[i];

      // Compute problem functions at trial point
      casadi_int info = evaluate(m, m->trial_xk, &objTrial, m->gk);
      m->nFunCalls++;
      cNormTrial = lInfConstraintNorm(m, m->trial_xk, m->gk);
      // Reduce step if evaluation fails, if lower bound is violated
      // or if objective or a constraint is NaN
      if (info != 0 || objTrial < obj_lo_ || objTrial > obj_up_
        || !(objTrial == objTrial) || !(cNormTrial == cNormTrial)) {
        print("info=%i, objTrial=%g\n", info, objTrial);
        // evaluation error, reduce stepsize
        reduceStepsize(m, &alpha);
        continue;
      } else {
        acceptStep(m, alpha);
        return 0;
      }
    }
    return 1;
  }

  /**
   *
   * Backtracking line search based on a filter
   * as described in Ipopt paper (Waechter 2006)
   *
   */
  casadi_int Blocksqp::filterLineSearch(BlocksqpMemory* m) const {
    auto d_nlp = &m->d_nlp;
    double alpha = 1.0;
    double cNormTrial=0, objTrial, dfTdeltaXi=0;

    // Compute ||constr(xi)|| at old point
    double cNorm = lInfConstraintNorm(m, d_nlp->z, m->gk);

    // Backtracking line search
    casadi_int k;
    for (k=0; k<max_line_search_; k++) {
      // Compute new trial point
      for (casadi_int i=0; i<nx_; i++)
        m->trial_xk[i] = d_nlp->z[i] + alpha * m->dxk[i];

      // Compute grad(f)^T * deltaXi
      dfTdeltaXi = 0.0;
      for (casadi_int i=0; i<nx_; i++)
        dfTdeltaXi += m->grad_fk[i] * m->dxk[i];

      // Compute objective and at ||constr(trial_xk)||_1 at trial point
      casadi_int info = evaluate(m, m->trial_xk, &objTrial, m->gk);
      m->nFunCalls++;
      cNormTrial = lInfConstraintNorm(m, m->trial_xk, m->gk);
      // Reduce step if evaluation fails, if lower bound is violated or if objective is NaN
      if (info != 0 || objTrial < obj_lo_ || objTrial > obj_up_
        || !(objTrial == objTrial) || !(cNormTrial == cNormTrial)) {
          // evaluation error, reduce stepsize
          reduceStepsize(m, &alpha);
          continue;
        }

      // Check acceptability to the filter
      if (pairInFilter(m, cNormTrial, objTrial)) {
        // Trial point is in the prohibited region defined by
        // the filter, try second order correction
        if (secondOrderCorrection(m, cNorm, cNormTrial, dfTdeltaXi, 0, k)) {
          break; // SOC yielded suitable alpha, stop
        } else {
          reduceStepsize(m, &alpha);
          continue;
        }
      }

      // Check sufficient decrease, case I:
      // If we are (almost) feasible and a "switching condition" is satisfied
      // require sufficient progress in the objective instead of bi-objective condition
      if (cNorm <= theta_min_) {
        // Switching condition, part 1: grad(f)^T * deltaXi < 0 ?
        if (dfTdeltaXi < 0)
          // Switching condition, part 2: alpha * (- grad(f)^T * deltaXi)**sF
          // > delta * cNorm**sTheta ?
          if (alpha * pow((-dfTdeltaXi), s_f_)
          > delta_ * pow(cNorm, s_theta_)) {
            // Switching conditions hold: Require satisfaction of Armijo condition for objective
            if (objTrial > m->obj + eta_*alpha*dfTdeltaXi) {
              // Armijo condition violated, try second order correction
              if (secondOrderCorrection(m, cNorm, cNormTrial, dfTdeltaXi, 1, k)) {
                break; // SOC yielded suitable alpha, stop
              } else {
                reduceStepsize(m, &alpha);
                continue;
              }
            } else {
              // found suitable alpha, stop
              acceptStep(m, alpha);
              break;
            }
          }
      }

      // Check sufficient decrease, case II:
      // Bi-objective (filter) condition
      if (cNormTrial < (1.0 - gamma_theta_) * cNorm
      || objTrial < m->obj - gamma_f_ * cNorm) {
        // found suitable alpha, stop
        acceptStep(m, alpha);
        break;
      } else {
        // Trial point is dominated by current point, try second order correction
        if (secondOrderCorrection(m, cNorm, cNormTrial, dfTdeltaXi, 0, k)) {
          break; // SOC yielded suitable alpha, stop
        } else {
          reduceStepsize(m, &alpha);
          continue;
        }
      }
    }

    // No step could be found by the line search
    if (k == max_line_search_) return 1;

    // Augment the filter if switching condition or Armijo condition does not hold
    if (dfTdeltaXi >= 0) {
      augmentFilter(m, cNormTrial, objTrial);
    } else if (alpha * pow((-dfTdeltaXi), s_f_) > delta_ * pow(cNorm, s_theta_)) {
      // careful with neg. exponents!
      augmentFilter(m, cNormTrial, objTrial);
    } else if (objTrial <= m->obj + eta_*alpha*dfTdeltaXi) {
      augmentFilter(m, cNormTrial, objTrial);
    }

    return 0;
  }


  /**
   *
   * Perform a second order correction step, i.e. solve the QP:
   *
   * min_d d^TBd + d^T grad_fk
   * s.t.  bl <= A^Td + constr(xi+alpha*deltaXi) - A^TdeltaXi <= bu
   *
   */
  bool Blocksqp::
  secondOrderCorrection(BlocksqpMemory* m, double cNorm, double cNormTrial,
    double dfTdeltaXi, bool swCond, casadi_int it) const {
    auto d_nlp = &m->d_nlp;

    // Perform SOC only on the first iteration of backtracking line search
    if (it > 0) return false;
    // If constraint violation of the trialstep is lower than the current one skip SOC
    if (cNormTrial < cNorm) return false;

    casadi_int nSOCS = 0;
    double cNormTrialSOC, cNormOld, objTrialSOC;

    // m->gk contains result at first trial point: c(xi+deltaXi)
    // m->jac_times_dxk and m->grad_fk are unchanged so far.

    // First SOC step
    std::vector<double> deltaXiSOC(nx_, 0.);
    std::vector<double> lambdaQPSOC(nx_+ng_, 0.);

    // Second order correction loop
    cNormOld = cNorm;
    for (casadi_int k=0; k<max_soc_iter_; k++) {
      nSOCS++;

      // Update bounds for SOC QP
      updateStepBounds(m, 1);

      // Solve SOC QP to obtain new, corrected deltaXi
      // (store in separate vector to avoid conflict with original deltaXi
      // -> need it in linesearch!)
      casadi_int info = solveQP(m, get_ptr(deltaXiSOC), get_ptr(lambdaQPSOC), false);
      if (info != 0) return false; // Could not solve QP, abort SOC

      // Set new SOC trial point
      for (casadi_int i=0; i<nx_; i++) {
        m->trial_xk[i] = d_nlp->z[i] + deltaXiSOC[i];
      }

      // Compute objective and ||constr(trialXiSOC)||_1 at SOC trial point
      info = evaluate(m, m->trial_xk, &objTrialSOC, m->gk);
      m->nFunCalls++;
      cNormTrialSOC = lInfConstraintNorm(m, m->trial_xk, m->gk);
      if (info != 0 || objTrialSOC < obj_lo_ || objTrialSOC > obj_up_
        || !(objTrialSOC == objTrialSOC) || !(cNormTrialSOC == cNormTrialSOC)) {
        return false; // evaluation error, abort SOC
      }

      // Check acceptability to the filter (in SOC)
      if (pairInFilter(m, cNormTrialSOC, objTrialSOC)) {
        // Trial point is in the prohibited region defined by the filter, abort SOC
        return false;
      }

      // Check sufficient decrease, case I (in SOC)
      // (Almost feasible and switching condition holds for line search alpha)
      if (cNorm <= theta_min_ && swCond) {
        if (objTrialSOC > m->obj + eta_*dfTdeltaXi) {
          // Armijo condition does not hold for SOC step, next SOC step

          // If constraint violation gets worse by SOC, abort
          if (cNormTrialSOC > kappa_soc_ * cNormOld) {
            return false;
          } else {
            cNormOld = cNormTrialSOC;
          }
          continue;
        } else {
          // found suitable alpha during SOC, stop
          acceptStep(m, get_ptr(deltaXiSOC), get_ptr(lambdaQPSOC), 1.0, nSOCS);
          return true;
        }
      }

      // Check sufficient decrease, case II (in SOC)
      if (cNorm > theta_min_ || !swCond) {
        if (cNormTrialSOC < (1.0 - gamma_theta_) * cNorm
        || objTrialSOC < m->obj - gamma_f_ * cNorm) {
          // found suitable alpha during SOC, stop
          acceptStep(m, get_ptr(deltaXiSOC), get_ptr(lambdaQPSOC), 1.0, nSOCS);
          return true;
        } else {
          // Trial point is dominated by current point, next SOC step

          // If constraint violation gets worse by SOC, abort
          if (cNormTrialSOC > kappa_soc_ * cNormOld) {
            return false;
          } else {
            cNormOld = cNormTrialSOC;
          }
          continue;
        }
      }
    }

    return false;
  }

  /**
   * Minimize constraint violation by solving an NLP with minimum norm objective
   *
   * "The dreaded restoration phase" -- Nick Gould
   */
  casadi_int Blocksqp::feasibilityRestorationPhase(BlocksqpMemory* m) const {
    auto d_nlp = &m->d_nlp;
    // No Feasibility restoration phase
    if (!restore_feas_) return -1;

    m->nRestPhaseCalls++;

    casadi_int ret, info;
    casadi_int maxRestIt = 100;
    casadi_int warmStart;
    double cNormTrial, objTrial, lStpNorm;


    // Create Input for the minimum norm NLP
    DMDict solver_in;

    // The reference point is the starting value for the restoration phase
    vector<double> in_x0(d_nlp->z, d_nlp->z+nx_);

    // Initialize slack variables such that the constraints are feasible
    for (casadi_int i=0; i<ng_; i++) {
      if (m->gk[i] <= d_nlp->lbz[i+nx_])
        in_x0.push_back(m->gk[i] - d_nlp->lbz[i+nx_]);
      else if (m->gk[i] > d_nlp->ubz[i+nx_])
        in_x0.push_back(m->gk[i] - d_nlp->ubz[i+nx_]);
      else
        in_x0.push_back(0.0);
    }

    // Add current iterate xk to parameter p
    vector<double> in_p(d_nlp->p, d_nlp->p+np_);
    vector<double> in_p2(d_nlp->z, d_nlp->z+nx_);
    in_p.insert(in_p.end(), in_p2.begin(), in_p2.end());

    // Set bounds for variables
    vector<double> in_lbx(d_nlp->lbz, d_nlp->lbz+nx_);
    vector<double> in_ubx(d_nlp->ubz, d_nlp->ubz+nx_);
    for (casadi_int i=0; i<ng_; i++) {
      in_lbx.push_back(-inf);
      in_ubx.push_back(inf);
    }

    // Set bounds for constraints
    vector<double> in_lbg(d_nlp->lbz+nx_, d_nlp->lbz+nx_+ng_);
    vector<double> in_ubg(d_nlp->ubz+nx_, d_nlp->ubz+nx_+ng_);

    solver_in["x0"] = in_x0;
    solver_in["p"] = in_p;
    solver_in["lbx"] = in_lbx;
    solver_in["ubx"] = in_ubx;
    solver_in["lbg"] = in_lbg;
    solver_in["ubg"] = in_ubg;

      /*

        Consider the following simple call:
        auto solver_out = rp_solver_(solver_in);

        This call in fact allocates memory,
        calls a memory-less eval(),
        and clears up the memory again.

        Since we want to access the memory later on,
        we need to unravel the simple call into its parts,
        and avoid the memory cleanup

      */

    // perform first iteration for the minimum norm NLP

    // Get the number of inputs and outputs
    int n_in = rp_solver_.n_in();
    int n_out = rp_solver_.n_out();

    // Get default inputs
    vector<DM> arg_v(n_in);
    for (casadi_int i=0; i<arg_v.size(); i++)
      arg_v[i] = DM::repmat(rp_solver_.default_in(i), rp_solver_.size1_in(i), 1);

    // Assign provided inputs
    for (auto&& e : solver_in) arg_v.at(rp_solver_.index_in(e.first)) = e.second;

    // Check sparsities
    for (casadi_int i=0; i<arg_v.size(); i++)
      casadi_assert_dev(arg_v[i].sparsity()==rp_solver_.sparsity_in(i));

    // Allocate results
    std::vector<DM> res(n_out);
    for (casadi_int i=0; i<n_out; i++) {
      if (res[i].sparsity()!=rp_solver_.sparsity_out(i))
        res[i] = DM::zeros(rp_solver_.sparsity_out(i));
    }

    // Allocate temporary memory if needed
    std::vector<casadi_int> iw_tmp(rp_solver_.sz_iw());
    std::vector<double> w_tmp(rp_solver_.sz_w());

    // Get pointers to input arguments
    std::vector<const double*> argp(rp_solver_.sz_arg());
    for (casadi_int i=0; i<n_in; ++i) argp[i] = get_ptr(arg_v[i]);

    // Get pointers to output arguments
    std::vector<double*> resp(rp_solver_.sz_res());
    for (casadi_int i=0; i<n_out; i++) resp[i] = get_ptr(res[i]);

    void* mem2 = rp_solver_.memory(0);

    // perform The m
    rp_solver_->eval(get_ptr(argp), get_ptr(resp), get_ptr(iw_tmp), get_ptr(w_tmp), mem2);

    // Get BlocksqpMemory and Blocksqp from restoration phase
    auto m2 = static_cast<BlocksqpMemory*>(mem2);
    const Blocksqp* bp = static_cast<const Blocksqp*>(rp_solver_.get());
    ret = m2->ret_;

    warmStart = 1;
    for (casadi_int it=0; it<maxRestIt; it++) {
      // One iteration for minimum norm NLP
      if (it > 0)
        ret = bp->run(m2, 1, warmStart);

      // If restMethod yields error, stop restoration phase
      if (ret == -1)
        break;

      // Get new xi from the restoration phase
      for (casadi_int i=0; i<nx_; i++)
        m->trial_xk[i] = m2->d_nlp.z[i];

      // Compute objective at trial point
      info = evaluate(m, m->trial_xk, &objTrial, m->gk);
      m->nFunCalls++;
      cNormTrial = lInfConstraintNorm(m, m->trial_xk, m->gk);
      if ( info != 0 || objTrial < obj_lo_ || objTrial > obj_up_ ||
          !(objTrial == objTrial) || !(cNormTrial == cNormTrial) )
        continue;

      // Is this iterate acceptable for the filter?
      if (!pairInFilter(m, cNormTrial, objTrial)) {
        // success
        print("Found a point acceptable for the filter.\n");
        ret = 0;
        break;
      }

      // If minimum norm NLP has converged, declare local infeasibility
        if (m2->tol < opttol_ && m2->cNormS < nlinfeastol_) {
          ret = 1;
          break;
        }
    }

    // Success or locally infeasible
    if (ret == 0 || ret == 1) {
        // Store the infinity norm of the multiplier step
        m->lambdaStepNorm = 0.0;
        // Compute restoration step
        for (casadi_int k=0; k<nx_; k++) {
            m->dxk[k] = d_nlp->z[k];

            d_nlp->z[k] = m->trial_xk[k];

            // Store lInf norm of dual step
            if ((lStpNorm = fabs(m2->lam_xk[k] - m->lam_xk[k])) > m->lambdaStepNorm)
                m->lambdaStepNorm = lStpNorm;
            m->lam_xk[k] = m2->lam_xk[k];
            m->lam_qp[k] = m2->lam_qp[k];

            m->dxk[k] -= d_nlp->z[k];
        }
        for (casadi_int k=0; k<ng_; k++) {
            // skip the dual variables for the slack variables in the restoration problem
            if ((lStpNorm = fabs(m2->lam_gk[k] - m->lam_gk[k])) > m->lambdaStepNorm)
                m->lambdaStepNorm = lStpNorm;
            m->lam_gk[k] = m2->lam_gk[k];
            m->lam_qp[k] = m2->lam_qp[nx_+ng_+k];
        }
        m->alpha = 1.0;
        m->nSOCS = 0;

        // reset reduced step counter
        m->reducedStepCount = 0;

        // reset Hessian and limited memory information
        resetHessian(m);
    }

    if (ret == 1) {
        if (print_iteration_) printProgress(m);
        updateStats(m);
        print("The problem seems to be locally infeasible. Infeasibilities minimized.\n");
    }

    return ret;
  }


  /**
   * Try to (partly) improve constraint violation by satisfying
   * the (pseudo) continuity constraints, i.e. do a single shooting
   * iteration with the current controls and measurement weights q and w
   */
  casadi_int Blocksqp::feasibilityRestorationHeuristic(BlocksqpMemory* m) const {
    auto d_nlp = &m->d_nlp;
    m->nRestHeurCalls++;

    // Call problem specific heuristic to reduce constraint violation.
    // For shooting methods that means setting consistent values for
    // shooting nodes by one forward integration.
    for (casadi_int k=0; k<nx_; k++) // input: last successful step
      m->trial_xk[k] = d_nlp->z[k];

    // FIXME(@jaeandersson) Not implemented
    return -1;
  }


  /**
   * If the line search fails, check if the full step reduces the KKT error by a factor kappaF.
   */
  casadi_int Blocksqp::kktErrorReduction(BlocksqpMemory* m) const {
    auto d_nlp = &m->d_nlp;
    casadi_int info = 0;
    double objTrial, cNormTrial, trialGradNorm, trialTol;

    // Compute new trial point
    for (casadi_int i=0; i<nx_; i++)
      m->trial_xk[i] = d_nlp->z[i] + m->dxk[i];

    // Compute objective and ||constr(trial_xk)|| at trial point
    std::vector<double> trialConstr(ng_, 0.);
    info = evaluate(m, m->trial_xk, &objTrial, get_ptr(trialConstr));
    m->nFunCalls++;
    cNormTrial = lInfConstraintNorm(m, m->trial_xk, get_ptr(trialConstr));
    if (info != 0 || objTrial < obj_lo_ || objTrial > obj_up_
      || !(objTrial == objTrial) || !(cNormTrial == cNormTrial)) {
      // evaluation error
      return 1;
    }

    // Compute KKT error of the new point

    // scaled norm of Lagrangian gradient
    std::vector<double> trialGradLagrange(nx_, 0.);
    calcLagrangeGradient(m, m->lam_qp, m->lam_qp+nx_, m->grad_fk,
                         m->jac_g,
                         get_ptr(trialGradLagrange), 0);

    trialGradNorm = casadi_norm_inf(nx_, get_ptr(trialGradLagrange));
    trialTol = trialGradNorm/(1.0+casadi_norm_inf(nx_+ng_, m->lam_qp));

    if (fmax(cNormTrial, trialTol) < kappa_f_ * fmax(m->cNorm, m->tol)) {
      acceptStep(m, 1.0);
      return 0;
    } else {
      return 1;
    }
  }

  /**
   * Check if current entry is accepted to the filter:
   * (cNorm, obj) in F_k
   */
  bool Blocksqp::
  pairInFilter(BlocksqpMemory* m, double cNorm, double obj) const {
    /*
     * A pair is in the filter if:
     * - it increases the objective and
     * - it also increases the constraint violation
     * The second expression in the if-clause states that we exclude
     * entries that are within the feasibility tolerance, e.g.
     * if an entry improves the constraint violation from 1e-16 to 1e-17,
     * but increases the objective considerably we also think of this entry
     * as dominated
     */

    for (auto&& f : m->filter) {
      if ((cNorm >= (1.0 - gamma_theta_) * f.first ||
           (cNorm < 0.01 * nlinfeastol_ && f.first < 0.01 * nlinfeastol_)) &&
          obj >= f.second - gamma_f_ * f.first) {
        return 1;
      }
    }

    return 0;
  }


  void Blocksqp::initializeFilter(BlocksqpMemory* m) const {
    std::pair<double, double> initPair(theta_max_, obj_lo_);

    // Remove all elements
    auto iter=m->filter.begin();
    while (iter != m->filter.end()) {
      std::set< std::pair<double, double> >::iterator iterToRemove = iter;
      iter++;
      m->filter.erase(iterToRemove);
    }

    // Initialize with pair (maxConstrViolation, objLowerBound);
    m->filter.insert(initPair);
  }

  /**
   * Augment the filter:
   * F_k+1 = F_k U { (c,f) | c > (1-gammaTheta)cNorm and f > obj-gammaF*c
   */
  void Blocksqp::
  augmentFilter(BlocksqpMemory* m, double cNorm, double obj) const {
    std::pair<double, double> entry((1-gamma_theta_)*cNorm,
                                    obj-gamma_f_*cNorm);

    // Augment filter by current element
    m->filter.insert(entry);

    // Remove dominated elements
    auto iter=m->filter.begin();
    while (iter != m->filter.end()) {
      if (iter->first > entry.first && iter->second > entry.second) {
        auto iterToRemove = iter;
        iter++;
        m->filter.erase(iterToRemove);
      } else {
        iter++;
      }
    }
  }

  /**
   * Initial Hessian: Identity matrix
   */
  void Blocksqp::calcInitialHessian(BlocksqpMemory* m) const {
    for (casadi_int b=0; b<nblocks_; b++)
      //if objective derv is computed exactly, don't set the last block!
      if (!(which_second_derv_ == 1 && block_hess_
        && b == nblocks_-1))
        calcInitialHessian(m, b);
  }


  /**
   * Initial Hessian for one block: Identity matrix
   */
  void Blocksqp::calcInitialHessian(BlocksqpMemory* m, casadi_int b) const {
    casadi_int dim = dim_[b];
    casadi_fill(m->hess[b], dim*dim, 0.);

    // Each block is a diagonal matrix
    for (casadi_int i=0; i<dim; i++)
      m->hess[b][i+i*dim] = ini_hess_diag_;

    // If we maintain 2 Hessians, also reset the second one
    if (m->hess2 != nullptr) {
      casadi_fill(m->hess2[b], dim*dim, 0.);
      for (casadi_int i=0; i<dim; i++)
        m->hess2[b][i+i*dim] = ini_hess_diag_;
    }
  }


  void Blocksqp::resetHessian(BlocksqpMemory* m) const {
    for (casadi_int b=0; b<nblocks_; b++) {
      if (!(which_second_derv_ == 1 && block_hess_ && b == nblocks_ - 1)) {
        // if objective derv is computed exactly, don't set the last block!
        resetHessian(m, b);
      }
    }
  }


  void Blocksqp::resetHessian(BlocksqpMemory* m, casadi_int b) const {
    casadi_int dim = dim_[b];

    // smallGamma and smallDelta are either subvectors of gamma and delta
    // or submatrices of gammaMat, deltaMat, i.e. subvectors of gamma and delta
    // from m prev. iterations (for L-BFGS)
    double *smallGamma = m->gammaMat + blocks_[b];
    double *smallDelta = m->deltaMat + blocks_[b];

    for (casadi_int i=0; i<hess_memsize_; ++i) {
      // Remove past information on Lagrangian gradient difference
      casadi_fill(smallGamma, dim, 0.);
      smallGamma += nx_;

      // Remove past information on steps
      smallDelta += nx_;
      casadi_fill(smallDelta, dim, 0.);
    }

    // Remove information on old scalars (used for COL sizing)
    m->delta_norm[b] = 1.0;
    m->delta_gamma[b] = 0.0;
    m->delta_norm_old[b] = 1.0;
    m->delta_gamma_old[b] = 0.0;

    m->noUpdateCounter[b] = -1;

    calcInitialHessian(m, b);
  }

  void Blocksqp::
  sizeInitialHessian(BlocksqpMemory* m, const double* gamma,
                     const double* delta, casadi_int b, casadi_int option) const {
    casadi_int dim = dim_[b];
    double scale;
    double myEps = 1.0e3 * eps_;

    if (option == 1) {
      // Shanno-Phua
      scale = casadi_dot(dim, gamma, gamma)
        / fmax(casadi_dot(dim, delta, gamma), myEps);
    } else if (option == 2) {
      // Oren-Luenberger
      scale = casadi_dot(dim, delta, gamma)
        / fmax(casadi_dot(dim, delta, delta), myEps);
      scale = fmin(scale, 1.0);
    } else if (option == 3) {
      // Geometric mean of 1 and 2
      scale = sqrt(casadi_dot(dim, gamma, gamma)
        / fmax(casadi_dot(dim, delta, delta), myEps));
    } else {
      // Invalid option, ignore
      return;
    }

    if (scale > 0.0) {
      scale = fmax(scale, myEps);
      for (casadi_int i=0; i<dim; i++)
        for (casadi_int j=0; j<dim; j++)
          m->hess[b][i+j*dim] *= scale;
    } else {
      scale = 1.0;
    }

    // statistics: average sizing factor
    m->averageSizingFactor += scale;
  }


  void Blocksqp::
  sizeHessianCOL(BlocksqpMemory* m, const double* gamma,
                 const double* delta, casadi_int b) const {
    casadi_int dim = dim_[b];
    double theta, scale, myEps = 1.0e3 * eps_;
    double deltaNorm, deltaNormOld, deltaGamma, deltaGammaOld, deltaBdelta;

    // Get sTs, sTs_, sTy, sTy_, sTBs
    deltaNorm = m->delta_norm[b];
    deltaGamma = m->delta_gamma[b];
    deltaNormOld = m->delta_norm_old[b];
    deltaGammaOld = m->delta_gamma_old[b];
    deltaBdelta = 0.0;
    for (casadi_int i=0; i<dim; i++)
      for (casadi_int j=0; j<dim; j++)
        deltaBdelta += delta[i] * m->hess[b][i+j*dim] * delta[j];

    // Centered Oren-Luenberger factor
    if (m->noUpdateCounter[b] == -1) {
      // in the first iteration, this should equal the OL factor
      theta = 1.0;
    } else {
      theta = fmin(col_tau1_, col_tau2_ * deltaNorm);
    }
    if (deltaNorm > myEps && deltaNormOld > myEps) {
      scale = (1.0 - theta)*deltaGammaOld / deltaNormOld + theta*deltaBdelta / deltaNorm;
      if (scale > eps_)
        scale = ((1.0 - theta)*deltaGammaOld / deltaNormOld + theta*deltaGamma / deltaNorm) / scale;
    } else {
      scale = 1.0;
    }

    // Size only if factor is between zero and one
    if (scale < 1.0 && scale > 0.0) {
      scale = fmax(col_eps_, scale);
      //print("Sizing value (COL) block %i = %g\n", b, scale);
      for (casadi_int i=0; i<dim; i++)
        for (casadi_int j=0; j<dim; j++)
          m->hess[b][i+j*dim] *= scale;

      // statistics: average sizing factor
      m->averageSizingFactor += scale;
    } else {
      m->averageSizingFactor += 1.0;
    }
  }

  /**
   * Apply BFGS or SR1 update blockwise and size blocks
   */
  void Blocksqp::
  calcHessianUpdate(BlocksqpMemory* m, casadi_int updateType, casadi_int hessScaling) const {
    casadi_int nBlocks;
    bool firstIter;

    //if objective derv is computed exactly, don't set the last block!
    if (which_second_derv_ == 1 && block_hess_)
      nBlocks = nblocks_ - 1;
    else
      nBlocks = nblocks_;

    // Statistics: how often is damping active, what is the average COL sizing factor?
    m->hessDamped = 0;
    m->averageSizingFactor = 0.0;

    for (casadi_int b=0; b<nBlocks; b++) {
      casadi_int dim = dim_[b];

      // smallGamma and smallDelta are subvectors of gamma and delta,
      // corresponding to partially separability
      double* smallGamma = m->gammaMat + blocks_[b];
      double* smallDelta = m->deltaMat + blocks_[b];

      // Is this the first iteration or the first after a Hessian reset?
      firstIter = (m->noUpdateCounter[b] == -1);

      // Update sTs, sTs_ and sTy, sTy_
      m->delta_norm_old[b] = m->delta_norm[b];
      m->delta_gamma_old[b] = m->delta_gamma[b];
      m->delta_norm[b] = casadi_dot(dim, smallDelta, smallDelta);
      m->delta_gamma[b] = casadi_dot(dim, smallDelta, smallGamma);

      // Sizing before the update
      if (hessScaling < 4 && firstIter)
        sizeInitialHessian(m, smallGamma, smallDelta, b, hessScaling);
      else if (hessScaling == 4)
        sizeHessianCOL(m, smallGamma, smallDelta, b);

      // Compute the new update
      if (updateType == 1) {
        calcSR1(m, smallGamma, smallDelta, b);

        // Prepare to compute fallback update as well
        m->hess = m->hess2;

        // Sizing the fallback update
        if (fallback_scaling_ < 4 && firstIter)
          sizeInitialHessian(m, smallGamma, smallDelta, b, fallback_scaling_);
        else if (fallback_scaling_ == 4)
          sizeHessianCOL(m, smallGamma, smallDelta, b);

        // Compute fallback update
        if (fallback_update_ == 2)
          calcBFGS(m, smallGamma, smallDelta, b);

        // Reset pointer
        m->hess = m->hess1;
      } else if (updateType == 2) {
        calcBFGS(m, smallGamma, smallDelta, b);
      }

      // If an update is skipped to often, reset Hessian block
      if (m->noUpdateCounter[b] > max_consec_skipped_updates_) {
        resetHessian(m, b);
      }
    }

    // statistics: average sizing factor
    m->averageSizingFactor /= nBlocks;
  }


  void Blocksqp::
  calcHessianUpdateLimitedMemory(BlocksqpMemory* m,
      casadi_int updateType, casadi_int hessScaling) const {
    casadi_int nBlocks;
    casadi_int m2, pos, posOldest, posNewest;
    casadi_int hessDamped, hessSkipped;
    double averageSizingFactor;

    //if objective derv is computed exactly, don't set the last block!
    if (which_second_derv_ == 1 && block_hess_) {
      nBlocks = nblocks_ - 1;
    } else {
      nBlocks = nblocks_;
    }

    // Statistics: how often is damping active, what is the average COL sizing factor?
    m->hessDamped = 0;
    m->hessSkipped = 0;
    m->averageSizingFactor = 0.0;

    for (casadi_int b=0; b<nBlocks; b++) {
      casadi_int dim = dim_[b];

      // smallGamma and smallDelta are submatrices of gammaMat, deltaMat,
      // i.e. subvectors of gamma and delta from m prev. iterations
      double *smallGamma = m->gammaMat + blocks_[b];
      double *smallDelta = m->deltaMat + blocks_[b];

      // Memory structure
      if (m->itCount > hess_memsize_) {
        m2 = hess_memsize_;
        posOldest = m->itCount % m2;
        posNewest = (m->itCount-1) % m2;
      } else {
        m2 = m->itCount;
        posOldest = 0;
        posNewest = m2-1;
      }

      // Set B_0 (pretend it's the first step)
      calcInitialHessian(m, b);
      m->delta_norm[b] = 1.0;
      m->delta_norm_old[b] = 1.0;
      m->delta_gamma[b] = 0.0;
      m->delta_gamma_old[b] = 0.0;
      m->noUpdateCounter[b] = -1;

      // Size the initial update, but with the most recent delta/gamma-pair
      double *gammai = smallGamma + nx_*posNewest;
      double *deltai = smallDelta + nx_*posNewest;
      sizeInitialHessian(m, gammai, deltai, b, hessScaling);

      for (casadi_int i=0; i<m2; i++) {
        pos = (posOldest+i) % m2;

        // Get new vector from list
        gammai = smallGamma + nx_*pos;
        deltai = smallDelta + nx_*pos;

        // Update sTs, sTs_ and sTy, sTy_
        m->delta_norm_old[b] = m->delta_norm[b];
        m->delta_gamma_old[b] = m->delta_gamma[b];
        m->delta_norm[b] = casadi_dot(dim, deltai, deltai);
        m->delta_gamma[b] = casadi_dot(dim, gammai, deltai);

        // Save statistics, we want to record them only for the most recent update
        averageSizingFactor = m->averageSizingFactor;
        hessDamped = m->hessDamped;
        hessSkipped = m->hessSkipped;

        // Selective sizing before the update
        if (hessScaling == 4) sizeHessianCOL(m, gammai, deltai, b);

        // Compute the new update
        if (updateType == 1) {
          calcSR1(m, gammai, deltai, b);
        } else if (updateType == 2) {
          calcBFGS(m, gammai, deltai, b);
        }

        m->nTotalUpdates++;

        // Count damping statistics only for the most recent update
        if (pos != posNewest) {
          m->hessDamped = hessDamped;
          m->hessSkipped = hessSkipped;
          if (hessScaling == 4)
            m->averageSizingFactor = averageSizingFactor;
        }
      }

      // If an update is skipped to often, reset Hessian block
      if (m->noUpdateCounter[b] > max_consec_skipped_updates_) {
        resetHessian(m, b);
      }
    }
    //blocks
    m->averageSizingFactor /= nBlocks;
  }


  void Blocksqp::
  calcHessianUpdateExact(BlocksqpMemory* m) const {
    // compute exact hessian
    (void)evaluate(m, m->exact_hess_lag);

    // assign exact hessian to blocks
    const casadi_int* col = exact_hess_lag_sp_.colind();
    const casadi_int* row = exact_hess_lag_sp_.row();
    casadi_int s, dim;

    for (casadi_int k=0; k<nblocks_; k++) {
      s = blocks_[k];
      dim = dim_[k];
      for (casadi_int i=0; i<dim; i++)
        // Set diagonal to 0 (may have been 1 because of CalcInitialHessian)
        m->hess[k][i + i*dim] = 0.0;
      for (casadi_int j=0; j<dim; j++) {
        for (casadi_int i=col[j+s]; i<col[j+1+s]; i++) {
          m->hess[k][row[i]-row[col[s]] + j*dim] = m->exact_hess_lag[i];
          if (row[i]-row[col[s]] < j)
            m->hess[k][j + (row[i]-row[col[s]])*dim] = m->exact_hess_lag[i];
        }
      }
    }

    // Prepare to compute fallback update as well
    m->hess = m->hess2;
    if (fallback_update_ == 2 && !hess_lim_mem_)
        calcHessianUpdate(m, fallback_update_, fallback_scaling_);
    else if (fallback_update_ == 0)
        calcInitialHessian(m);  // Set fallback as Identity

    // Reset pointer
    m->hess = m->hess1;
  }


  void Blocksqp::
  calcBFGS(BlocksqpMemory* m, const double* gamma,
    const double* delta, casadi_int b) const {
    casadi_int dim = dim_[b];
    double h1 = 0.0;
    double h2 = 0.0;
    double thetaPowell = 0.0;
    casadi_int damped;

    /* Work with a local copy of gamma because damping may need to change gamma.
     * Note that m->gamma needs to remain unchanged!
     * This may be important in a limited memory context:
     * When information is "forgotten", B_i-1 is different and the
     *  original gamma might lead to an undamped update with the new B_i-1! */
    std::vector<double> gamma2(gamma, gamma+dim);

    double *B = m->hess[b];

    // Bdelta = B*delta (if sizing is enabled, B is the sized B!)
    // h1 = delta^T * B * delta
    // h2 = delta^T * gamma
    vector<double> Bdelta(dim, 0.0);
    for (casadi_int i=0; i<dim; i++) {
      for (casadi_int k=0; k<dim; k++)
        Bdelta[i] += B[i+k*dim] * delta[k];

      h1 += delta[i] * Bdelta[i];
      //h2 += delta[i] * gamma[i];
    }
    h2 = m->delta_gamma[b];

    /* Powell's damping strategy to maintain pos. def. (Nocedal/Wright p.537; SNOPT paper)
     * Interpolates between current approximation and unmodified BFGS */
    damped = 0;
    if (hess_damp_)
      if (h2 < hess_damp_fac_ * h1 / m->alpha && fabs(h1 - h2) > 1.0e-12) {
        // At the first iteration h1 and h2 are equal due to COL scaling

        thetaPowell = (1.0 - hess_damp_fac_)*h1 / (h1 - h2);

        // Redefine gamma and h2 = delta^T * gamma
        h2 = 0.0;
        for (casadi_int i=0; i<dim; i++) {
          gamma2[i] = thetaPowell*gamma2[i] + (1.0 - thetaPowell)*Bdelta[i];
          h2 += delta[i] * gamma2[i];
        }

        // Also redefine deltaGamma for computation of sizing factor in the next iteration
        m->delta_gamma[b] = h2;

        damped = 1;
      }

    // For statistics: count number of damped blocks
    m->hessDamped += damped;

    // B_k+1 = B_k - Bdelta * (Bdelta)^T / h1 + gamma * gamma^T / h2
    double myEps = 1.0e2 * eps_;
    if (fabs(h1) < myEps || fabs(h2) < myEps) {
      // don't perform update because of bad condition, might introduce negative eigenvalues
      m->noUpdateCounter[b]++;
      m->hessDamped -= damped;
      m->hessSkipped++;
      m->nTotalSkippedUpdates++;
    } else {
      for (casadi_int i=0; i<dim; i++)
        for (casadi_int j=0; j<dim; j++)
          B[i+j*dim] += - Bdelta[i]*Bdelta[j]/h1 + gamma2[i]*gamma2[j]/h2;

      m->noUpdateCounter[b] = 0;
    }
  }


  void Blocksqp::
  calcSR1(BlocksqpMemory* m, const double* gamma,
          const double* delta, casadi_int b) const {
    casadi_int dim = dim_[b];
    double *B = m->hess[b];
    double myEps = 1.0e2 * eps_;
    double r = 1.0e-8;
    double h = 0.0;

    // gmBdelta = gamma - B*delta
    // h = (gamma - B*delta)^T * delta
    vector<double> gmBdelta(dim);
    for (casadi_int i=0; i<dim; i++) {
      gmBdelta[i] = gamma[i];
      for (casadi_int k=0; k<dim; k++)
        gmBdelta[i] -= B[i+k*dim] * delta[k];

      h += (gmBdelta[i] * delta[i]);
    }

    // B_k+1 = B_k + gmBdelta * gmBdelta^T / h
    if (fabs(h) < r * casadi_norm_2(dim, delta)
      *casadi_norm_2(dim, get_ptr(gmBdelta)) || fabs(h) < myEps) {
      // Skip update if denominator is too small
      m->noUpdateCounter[b]++;
      m->hessSkipped++;
      m->nTotalSkippedUpdates++;
    } else {
      for (casadi_int i=0; i<dim; i++)
        for (casadi_int j=0; j<dim; j++)
          B[i+j*dim] += gmBdelta[i]*gmBdelta[j]/h;
      m->noUpdateCounter[b] = 0;
    }
  }


  /**
   * Set deltaXi and gamma as a column in the matrix containing
   * the m most recent delta and gamma
   */
  void Blocksqp::updateDeltaGamma(BlocksqpMemory* m) const {
    if (hess_memsize_ == 1) return;

    m->dxk = m->deltaMat + nx_*(m->itCount % hess_memsize_);
    m->gamma = m->gammaMat + nx_*(m->itCount % hess_memsize_);
  }

  void Blocksqp::
  computeNextHessian(BlocksqpMemory* m, casadi_int idx, casadi_int maxQP) const {
    // Compute fallback update only once
    if (idx == 1) {
        // Switch storage
        m->hess = m->hess2;

        // If last block contains exact Hessian, we need to copy it
        if (which_second_derv_ == 1) {
          casadi_int dim = dim_[nblocks_-1];
          casadi_copy(m->hess1[nblocks_-1], dim*dim, m->hess2[nblocks_-1]);
        }

        // Limited memory: compute fallback update only when needed
        if (hess_lim_mem_) {
            m->itCount--;
            casadi_int hessDampSave = hess_damp_;
            const_cast<Blocksqp*>(this)->hess_damp_ = 1;
            calcHessianUpdateLimitedMemory(m, fallback_update_, fallback_scaling_);
            const_cast<Blocksqp*>(this)->hess_damp_ = hessDampSave;
            m->itCount++;
          }
        /* Full memory: both updates must be computed in every iteration
         * so switching storage is enough */
      }

    // 'Nontrivial' convex combinations
    if (maxQP > 2) {
        /* Convexification parameter: mu_l = l / (maxQP-1).
         * Compute it only in the first iteration, afterwards update
         * by recursion: mu_l/mu_(l-1) */
        double idxF = idx;
        double mu = (idx==1) ? 1.0 / (maxQP-1) : idxF / (idxF - 1.0);
        double mu1 = 1.0 - mu;
        for (casadi_int b=0; b<nblocks_; b++) {
          casadi_int dim = dim_[b];
          for (casadi_int i=0; i<dim; i++) {
            for (casadi_int j=0; j<dim; j++) {
              m->hess2[b][i+j*dim] *= mu;
              m->hess2[b][i+j*dim] += mu1 * m->hess1[b][i+j*dim];
            }
          }
        }
    }
  }


  /**
   * Inner loop of SQP algorithm:
   * Solve a sequence of QPs until pos. def. assumption (G3*) is satisfied.
   */
  casadi_int Blocksqp::
  solveQP(BlocksqpMemory* m, double* deltaXi, double* lambdaQP,
    bool matricesChanged) const {
    casadi_int maxQP, l;
    if (globalization_ &&
        (hess_update_ == 1 || hess_update_ == 4) &&
        matricesChanged &&
        m->itCount > 1) {
        maxQP = max_conv_qp_ + 1;
      } else {
      maxQP = 1;
    }

    /*
     * Prepare for qpOASES
     */

    // Setup QProblem data
    if (matricesChanged) {
      delete m->A;
      m->A = nullptr;
      copy_vector(Asp_.colind(), m->colind);
      copy_vector(Asp_.row(), m->row);
      int* jacIndRow = get_ptr(m->row);
      int* jacIndCol = get_ptr(m->colind);
      m->A = new qpOASES::SparseMatrix(ng_, nx_,
                                       jacIndRow, jacIndCol, m->jac_g);
    }
    double *g = m->grad_fk;
    double *lb = m->lbx_qp;
    double *lu = m->ubx_qp;
    double *lbA = m->lba_qp;
    double *luA = m->uba_qp;

    // qpOASES options
    qpOASES::Options opts;
    if (matricesChanged && maxQP > 1)
      opts.enableInertiaCorrection = qpOASES::BT_FALSE;
    else
      opts.enableInertiaCorrection = qpOASES::BT_TRUE;
    opts.enableEqualities = qpOASES::BT_TRUE;
    opts.initialStatusBounds = qpOASES::ST_INACTIVE;
    opts.printLevel = qpOASES::PL_NONE;
    opts.numRefinementSteps = 2;
    opts.epsLITests =  2.2204e-08;
    m->qp->setOptions(opts);

    // Other variables for qpOASES
    double cpuTime = matricesChanged ? max_time_qp_ : 0.1*max_time_qp_;
    int maxIt = matricesChanged ? max_it_qp_ : static_cast<int>(0.1*max_it_qp_);
    qpOASES::SolutionAnalysis solAna;
    qpOASES::returnValue ret = qpOASES::RET_INFO_UNDEFINED;

    /*
     * QP solving loop for convex combinations (sequential)
     */
    for (l=0; l<maxQP; l++) {
        /*
         * Compute a new Hessian
         */
        if (l > 0) {
          // If the solution of the first QP was rejected, consider second Hessian
          m->qpResolve++;
          computeNextHessian(m, l, maxQP);
        }

        if (l == maxQP-1) {
          // Enable inertia correction for supposedly convex QPs, just in case
          opts.enableInertiaCorrection = qpOASES::BT_TRUE;
          m->qp->setOptions(opts);
        }

        /*
         * Prepare the current Hessian for qpOASES
         */
        if (matricesChanged) {
          // Convert block-Hessian to sparse format
          convertHessian(m);
          delete m->H;
          m->H = nullptr;
          m->H = new qpOASES::SymSparseMat(nx_, nx_,
                                           m->hessIndRow, m->hessIndCol,
                                           m->hess_lag);
          m->H->createDiagInfo();
        }

        /*
         * Call qpOASES
         */
        if (matricesChanged) {
            maxIt = max_it_qp_;
            cpuTime = max_time_qp_;
            if (m->qp->getStatus() == qpOASES::QPS_HOMOTOPYQPSOLVED ||
                m->qp->getStatus() == qpOASES::QPS_SOLVED) {
                ret = m->qp->hotstart(m->H, g, m->A, lb, lu, lbA, luA, maxIt, &cpuTime);
              } else {
                if (warmstart_) {
                  ret = m->qp->init(m->H, g, m->A, lb, lu, lbA, luA, maxIt, &cpuTime,
                    deltaXi, lambdaQP);
                } else {
                  ret = m->qp->init(m->H, g, m->A, lb, lu, lbA, luA, maxIt, &cpuTime);
                }
              }
          } else {
            // Second order correction: H and A do not change
            maxIt = static_cast<int>(0.1*max_it_qp_);
            cpuTime = 0.1*max_time_qp_;
            ret = m->qp->hotstart(g, lb, lu, lbA, luA, maxIt, &cpuTime);
          }

        /*
         * Check assumption (G3*) if nonconvex QP was solved
         */
        if (l < maxQP-1 && matricesChanged) {
            if (ret == qpOASES::SUCCESSFUL_RETURN) {
                if (schur_) {
                  ret = solAna.checkCurvatureOnStronglyActiveConstraints(
                    dynamic_cast<qpOASES::SQProblemSchur*>(m->qp));
                } else {
                  ret = solAna.checkCurvatureOnStronglyActiveConstraints(m->qp);
                }
              }

            if (ret == qpOASES::SUCCESSFUL_RETURN) {
              // QP was solved successfully and curvature is positive after removing bounds
                m->qpIterations = maxIt + 1;
                break; // Success!
              } else {
              // QP solution is rejected, save statistics
                if (ret == qpOASES::RET_SETUP_AUXILIARYQP_FAILED)
                  m->qpIterations2++;
                else
                  m->qpIterations2 += maxIt + 1;
                m->rejectedSR1++;
              }
          } else {
            // Convex QP was solved, no need to check assumption (G3*)
            m->qpIterations += maxIt + 1;
          }

      } // End of QP solving loop

    /*
     * Post-processing
     */

    // Get solution from qpOASES
    m->qp->getPrimalSolution(deltaXi);
    m->qp->getDualSolution(lambdaQP);
    m->qpObj = m->qp->getObjVal();

    // Compute constrJac*deltaXi, need this for second order correction step
    casadi_fill(m->jac_times_dxk, ng_, 0.);
    casadi_mv(m->jac_g, Asp_, deltaXi, m->jac_times_dxk, 0);

    // Print qpOASES error code, if any
    if (ret != qpOASES::SUCCESSFUL_RETURN && matricesChanged)
      print("***WARNING: qpOASES error message: \"%s\"\n",
              qpOASES::MessageHandling::getErrorCodeMessage(ret));

    // Point Hessian again to the first Hessian
    m->hess = m->hess1;

    /* For full-memory Hessian: Restore fallback Hessian if convex combinations
     * were used during the loop */
    if (!hess_lim_mem_ && maxQP > 2 && matricesChanged) {
        double mu = 1.0 / l;
        double mu1 = 1.0 - mu;
        casadi_int nBlocks = (which_second_derv_ == 1) ? nblocks_-1 : nblocks_;
        for (casadi_int b=0; b<nBlocks; b++) {
          casadi_int dim = dim_[b];
          for (casadi_int i=0; i<dim; i++) {
            for (casadi_int j=0; j<dim; j++) {
              m->hess2[b][i+j*dim] *= mu;
              m->hess2[b][i+j*dim] += mu1 * m->hess1[b][i+j*dim];
            }
          }
        }
    }

    /* Return code depending on qpOASES returnvalue
     * 0: Success
     * 1: Maximum number of iterations reached
     * 2: Unbounded
     * 3: Infeasible
     * 4: Other error */
     switch (ret) {
        case qpOASES::SUCCESSFUL_RETURN:
          return 0;
        case qpOASES::RET_MAX_NWSR_REACHED:
          return 1;
        case qpOASES::RET_HESSIAN_NOT_SPD:
        case qpOASES::RET_HESSIAN_INDEFINITE:
        case qpOASES::RET_INIT_FAILED_UNBOUNDEDNESS:
        case qpOASES::RET_QP_UNBOUNDED:
        case qpOASES::RET_HOTSTART_STOPPED_UNBOUNDEDNESS:
          return 2;
        case qpOASES::RET_INIT_FAILED_INFEASIBILITY:
        case qpOASES::RET_QP_INFEASIBLE:
        case qpOASES::RET_HOTSTART_STOPPED_INFEASIBILITY:
          return 3;
        default:
          return 4;
     }
  }

  /**
   * Set bounds on the step (in the QP), either according
   * to variable bounds in the NLP or according to
   * trust region box radius
   */
  void Blocksqp::updateStepBounds(BlocksqpMemory* m, bool soc) const {
    auto d_nlp = &m->d_nlp;
    // Bounds on step
    for (casadi_int i=0; i<nx_; i++) {
      double lbx = d_nlp->lbz[i];
      if (lbx != inf) {
        m->lbx_qp[i] = lbx - d_nlp->z[i];
      } else {
        m->lbx_qp[i] = inf;
      }

      double ubx = d_nlp->ubz[i];
      if (ubx != inf) {
        m->ubx_qp[i] = ubx - d_nlp->z[i];
      } else {
        m->ubx_qp[i] = inf;
      }
    }

    // Bounds on linearized constraints
    for (casadi_int i=0; i<ng_; i++) {
      double lbg = d_nlp->lbz[i+nx_];
      if (lbg != inf) {
        m->lba_qp[i] = lbg - m->gk[i];
        if (soc) m->lba_qp[i] += m->jac_times_dxk[i];
      } else {
        m->lba_qp[i] = inf;
      }

      double ubg = d_nlp->ubz[i+nx_];
      if (ubg != inf) {
        m->uba_qp[i] = ubg - m->gk[i];
        if (soc) m->uba_qp[i] += m->jac_times_dxk[i];
      } else {
        m->uba_qp[i] = inf;
      }
    }
  }

  void Blocksqp::printProgress(BlocksqpMemory* m) const {
    /*
     * m->steptype:
     *-1: full step was accepted because it reduces the KKT error although line search failed
     * 0: standard line search step
     * 1: Hessian has been reset to identity
     * 2: feasibility restoration heuristic has been called
     * 3: feasibility restoration phase has been called
     */

     // Print headline every twenty iterations
    if (m->itCount % 20 == 0) {
      print("%-8s", "   it");
      print("%-21s", " qpIt");
      print("%-9s", "obj");
      print("%-11s", "feas");
      print("%-7s", "opt");
      print("%-11s", "|lgrd|");
      print("%-9s", "|stp|");
      print("%-10s", "|lstp|");
      print("%-8s", "alpha");
      print("%-6s", "nSOCS");
      print("%-18s", "sk, da, sca");
      print("%-6s", "QPr,mu");
      print("\n");
    }

    if (m->itCount == 0) {
      // Values for first iteration
      print("%5i  ", m->itCount);
      print("%11i ", 0);
      print("% 10e  ", m->obj);
      print("%-10.2e", m->cNormS);
      print("%-10.2e", m->tol);
      print("\n");
    } else {
      // All values
      print("%5i  ", m->itCount);
      print("%5i+%5i ", m->qpIterations, m->qpIterations2);
      print("% 10e  ", m->obj);
      print("%-10.2e", m->cNormS);
      print("%-10.2e", m->tol);
      print("%-10.2e", m->gradNorm);
      print("%-10.2e", casadi_norm_inf(nx_, m->dxk));
      print("%-10.2e", m->lambdaStepNorm);
      print("%-9.1e", m->alpha);
      print("%5i", m->nSOCS);
      print("%3i, %3i, %-9.1e", m->hessSkipped, m->hessDamped, m->averageSizingFactor);
      print("%i, %-9.1e", m->qpResolve, casadi_norm_1(nblocks_, m->delta_h)/nblocks_);
      print("\n");
    }
  }

  void Blocksqp::initStats(BlocksqpMemory* m) const {
    m->itCount = 0;
    m->qpItTotal = 0;
    m->qpIterations = 0;
    m->hessSkipped = 0;
    m->hessDamped = 0;
    m->averageSizingFactor = 0.0;
  }

  void Blocksqp::updateStats(BlocksqpMemory* m) const {
    // Do not accidentally print hessSkipped in the next iteration
    m->hessSkipped = 0;
    m->hessDamped = 0;

    // qpIterations = number of iterations for the QP that determines the step,
    // can be a resolve (+SOC)
    // qpIterations2 = number of iterations for a QP which solution was discarded
    m->qpItTotal += m->qpIterations;
    m->qpItTotal += m->qpIterations2;
    m->qpIterations = 0;
    m->qpIterations2 = 0;
    m->qpResolve = 0;
  }

  /**
   * Allocate memory for variables
   * required by all optimization
   * algorithms except for the Jacobian
   */
  void Blocksqp::reset_sqp(BlocksqpMemory* m) const {
    // dual variables (for general constraints and variable bounds)
    casadi_fill(m->lam_xk, nx_, 0.);
    casadi_fill(m->lam_gk, ng_, 0.);

    // constraint vector with lower and upper bounds
    // (Box constraints are not included in the constraint list)
    casadi_fill(m->gk, ng_, 0.);

    // gradient of objective
    casadi_fill(m->grad_fk, nx_, 0.);

    // gradient of Lagrangian
    casadi_fill(m->grad_lagk, nx_, 0.);

    // current step
    casadi_fill(m->deltaMat, nx_*hess_memsize_, 0.);
    m->dxk = m->deltaMat;

    // trial step (temporary variable, for line search)
    casadi_fill(m->trial_xk, nx_, 0.);

    // bounds for step (QP subproblem)
    casadi_fill(m->lbx_qp, nx_, 0.);
    casadi_fill(m->ubx_qp, nx_, 0.);
    casadi_fill(m->lba_qp, ng_, 0.);
    casadi_fill(m->uba_qp, ng_, 0.);

    // product of constraint Jacobian with step (deltaXi)
    casadi_fill(m->jac_times_dxk, ng_, 0.);

    // dual variables of QP (simple bounds and general constraints)
    casadi_fill(m->lam_qp, nx_+ng_, 0.);

    // line search parameters
    casadi_fill(m->delta_h, nblocks_, 0.);

    // filter as a set of pairs
    m->filter.clear();

    // difference of Lagrangian gradients
    casadi_fill(m->gammaMat, nx_*hess_memsize_, 0.);
    m->gamma = m->gammaMat;

    // Scalars that are used in various Hessian update procedures
    casadi_fill(m->noUpdateCounter, nblocks_, casadi_int(-1));

    // For selective sizing: for each block save sTs, sTs_, sTy, sTy_
    casadi_fill(m->delta_norm, nblocks_, 1.);
    casadi_fill(m->delta_norm_old, nblocks_, 1.);
    casadi_fill(m->delta_gamma, nblocks_, 0.);
    casadi_fill(m->delta_gamma_old, nblocks_, 0.);

    // Create one Matrix for one diagonal block in the Hessian
    for (casadi_int b=0; b<nblocks_; b++) {
      casadi_int dim = dim_[b];
      casadi_fill(m->hess1[b], dim*dim, 0.);
    }

    // For SR1 or finite differences, maintain two Hessians
    if (hess_update_ == 1 || hess_update_ == 4) {
      for (casadi_int b=0; b<nblocks_; b++) {
        casadi_int dim = dim_[b];
        casadi_fill(m->hess2[b], dim*dim, 0.);
      }
    }

    // Set Hessian pointer
    m->hess = m->hess1;
  }

  /**
   * Convert array *hess to a single symmetric sparse matrix in
   * Harwell-Boeing format (as used by qpOASES)
   */
  void Blocksqp::
  convertHessian(BlocksqpMemory* m) const {
    casadi_int count, colCountTotal, rowOffset;
    casadi_int nnz;

    // 1) count nonzero elements
    nnz = 0;
    for (casadi_int b=0; b<nblocks_; b++) {
      casadi_int dim = dim_[b];
      for (casadi_int i=0; i<dim; i++) {
        for (casadi_int j=0; j<dim; j++) {
          if (fabs(m->hess[b][i+j*dim]) > eps_) {
            nnz++;
          }
        }
      }
    }

    m->hessIndCol = m->hessIndRow + nnz;
    m->hessIndLo = m->hessIndCol + (nx_+1);

    // 2) store matrix entries columnwise in hessNz
    count = 0; // runs over all nonzero elements
    colCountTotal = 0; // keep track of position in large matrix
    rowOffset = 0;
    for (casadi_int b=0; b<nblocks_; b++) {
      casadi_int dim = dim_[b];

      for (casadi_int i=0; i<dim; i++) {
        // column 'colCountTotal' starts at element 'count'
        m->hessIndCol[colCountTotal] = count;

        for (casadi_int j=0; j<dim; j++) {
          if (fabs(m->hess[b][i+j*dim]) > eps_) {
              m->hess_lag[count] = m->hess[b][i+j*dim];
              m->hessIndRow[count] = j + rowOffset;
              count++;
          }
        }
        colCountTotal++;
      }
      rowOffset += dim;
    }
    m->hessIndCol[colCountTotal] = count;

    // 3) Set reference to lower triangular matrix
    for (casadi_int j=0; j<nx_; j++) {
      casadi_int i;
      for (i=m->hessIndCol[j]; i<m->hessIndCol[j+1] && m->hessIndRow[i]<j; i++) {}
      m->hessIndLo[j] = i;
    }

    if (count != nnz)
      print("***WARNING: Error in convertHessian: %i elements processed, "
            "should be %i elements!\n", count, nnz);
  }

  void Blocksqp::initIterate(BlocksqpMemory* m) const {
    m->alpha = 1.0;
    m->nSOCS = 0;
    m->reducedStepCount = 0;
    m->steptype = 0;

    m->obj = inf;
    m->tol = inf;
    m->cNorm = theta_max_;
    m->gradNorm = inf;
    m->lambdaStepNorm = 0.0;
  }

  casadi_int Blocksqp::
  evaluate(BlocksqpMemory* m,
           double *f, double *g,
           double *grad_f, double *jac_g) const {
    auto d_nlp = &m->d_nlp;
    m->arg[0] = d_nlp->z; // x
    m->arg[1] = d_nlp->p; // p
    m->res[0] = f; // f
    m->res[1] = g; // g
    m->res[2] = grad_f; // grad:f:x
    m->res[3] = jac_g; // jac:g:x
    return calc_function(m, "nlp_gf_jg");
  }

  casadi_int Blocksqp::
  evaluate(BlocksqpMemory* m, const double *xk, double *f,
           double *g) const {
    auto d_nlp = &m->d_nlp;
    m->arg[0] = xk; // x
    m->arg[1] = d_nlp->p; // p
    m->res[0] = f; // f
    m->res[1] = g; // g
    return calc_function(m, "nlp_fg");
  }

  casadi_int Blocksqp::
  evaluate(BlocksqpMemory* m,
           double *exact_hess_lag) const {
    auto d_nlp = &m->d_nlp;
    static std::vector<double> ones;
    ones.resize(nx_);
    for (casadi_int i=0; i<nx_; ++i) ones[i] = 1.0;
    static std::vector<double> minus_lam_gk;
    minus_lam_gk.resize(ng_);
    // Langrange function in blocksqp is L = f - lambdaT * g, whereas + in casadi
    for (casadi_int i=0; i<ng_; ++i) minus_lam_gk[i] = -m->lam_gk[i];
    m->arg[0] = d_nlp->z; // x
    m->arg[1] = d_nlp->p; // p
    m->arg[2] = get_ptr(ones); // lam:f
    m->arg[3] = get_ptr(minus_lam_gk); // lam:g
    m->res[0] = exact_hess_lag; // hess:gamma:x:x
    return calc_function(m, "nlp_hess_l");;
  }

  BlocksqpMemory::BlocksqpMemory() {
    qpoases_mem = nullptr;
    H = nullptr;
    A = nullptr;
    qp = nullptr;
  }

  BlocksqpMemory::~BlocksqpMemory() {
    delete qpoases_mem;
    delete H;
    delete A;
    delete qp;
  }

  double Blocksqp::
  lInfConstraintNorm(BlocksqpMemory* m, const double* xk, const double* g) const {
    auto d_nlp = &m->d_nlp;
    return fmax(casadi_max_viol(nx_, xk, d_nlp->lbz, d_nlp->ubz),
                casadi_max_viol(ng_, g, d_nlp->lbz+nx_, d_nlp->ubz+nx_));
  }


  Blocksqp::Blocksqp(DeserializingStream& s) : Nlpsol(s) {
    s.version("Blocksqp", 1);
    s.unpack("Blocksqp::nblocks", nblocks_);
    s.unpack("Blocksqp::blocks", blocks_);
    s.unpack("Blocksqp::dim", dim_);
    s.unpack("Blocksqp::nnz_H", nnz_H_);
    s.unpack("Blocksqp::Asp", Asp_);
    s.unpack("Blocksqp::Hsp", Hsp_);
    s.unpack("Blocksqp::exact_hess_lag_sp_", exact_hess_lag_sp_);
    s.unpack("Blocksqp::linsol_plugin", linsol_plugin_);
    s.unpack("Blocksqp::print_header", print_header_);
    s.unpack("Blocksqp::print_iteration", print_iteration_);
    s.unpack("Blocksqp::eps", eps_);
    s.unpack("Blocksqp::opttol", opttol_);
    s.unpack("Blocksqp::nlinfeastol", nlinfeastol_);
    s.unpack("Blocksqp::schur", schur_);
    s.unpack("Blocksqp::globalization", globalization_);
    s.unpack("Blocksqp::restore_feas", restore_feas_);
    s.unpack("Blocksqp::max_line_search", max_line_search_);
    s.unpack("Blocksqp::max_consec_reduced_steps", max_consec_reduced_steps_);
    s.unpack("Blocksqp::max_consec_skipped_updates", max_consec_skipped_updates_);
    s.unpack("Blocksqp::max_it_qp", max_it_qp_);
    s.unpack("Blocksqp::max_iter", max_iter_);
    s.unpack("Blocksqp::warmstart", warmstart_);
    s.unpack("Blocksqp::qp_init", qp_init_);
    s.unpack("Blocksqp::block_hess", block_hess_);
    s.unpack("Blocksqp::hess_scaling", hess_scaling_);
    s.unpack("Blocksqp::fallback_scaling", fallback_scaling_);
    s.unpack("Blocksqp::max_time_qp", max_time_qp_);
    s.unpack("Blocksqp::ini_hess_diag", ini_hess_diag_);
    s.unpack("Blocksqp::col_eps", col_eps_);
    s.unpack("Blocksqp::col_tau1", col_tau1_);
    s.unpack("Blocksqp::col_tau2", col_tau2_);
    s.unpack("Blocksqp::hess_damp", hess_damp_);
    s.unpack("Blocksqp::hess_damp_fac", hess_damp_fac_);
    s.unpack("Blocksqp::hess_update", hess_update_);
    s.unpack("Blocksqp::fallback_update", fallback_update_);
    s.unpack("Blocksqp::hess_lim_mem", hess_lim_mem_);
    s.unpack("Blocksqp::hess_memsize", hess_memsize_);
    s.unpack("Blocksqp::which_second_derv", which_second_derv_);
    s.unpack("Blocksqp::skip_first_globalization", skip_first_globalization_);
    s.unpack("Blocksqp::conv_strategy", conv_strategy_);
    s.unpack("Blocksqp::max_conv_qp", max_conv_qp_);
    s.unpack("Blocksqp::max_soc_iter", max_soc_iter_);
    s.unpack("Blocksqp::gamma_theta", gamma_theta_);
    s.unpack("Blocksqp::gamma_f", gamma_f_);
    s.unpack("Blocksqp::kappa_soc", kappa_soc_);
    s.unpack("Blocksqp::kappa_f", kappa_f_);
    s.unpack("Blocksqp::theta_max", theta_max_);
    s.unpack("Blocksqp::theta_min", theta_min_);
    s.unpack("Blocksqp::delta", delta_);
    s.unpack("Blocksqp::s_theta", s_theta_);
    s.unpack("Blocksqp::s_f", s_f_);
    s.unpack("Blocksqp::kappa_minus", kappa_minus_);
    s.unpack("Blocksqp::kappa_plus", kappa_plus_);
    s.unpack("Blocksqp::kappa_plus_max", kappa_plus_max_);
    s.unpack("Blocksqp::delta_h0", delta_h0_);
    s.unpack("Blocksqp::eta", eta_);
    s.unpack("Blocksqp::obj_lo", obj_lo_);
    s.unpack("Blocksqp::obj_up", obj_up_);
    s.unpack("Blocksqp::rho", rho_);
    s.unpack("Blocksqp::zeta", zeta_);
    s.unpack("Blocksqp::rp_solver", rp_solver_);
    s.unpack("Blocksqp::print_maxit_reached", print_maxit_reached_);

  }

  void Blocksqp::serialize_body(SerializingStream &s) const {
    Nlpsol::serialize_body(s);
    s.version("Blocksqp", 1);
    s.pack("Blocksqp::nblocks", nblocks_);
    s.pack("Blocksqp::blocks", blocks_);
    s.pack("Blocksqp::dim", dim_);
    s.pack("Blocksqp::nnz_H", nnz_H_);
    s.pack("Blocksqp::Asp", Asp_);
    s.pack("Blocksqp::Hsp", Hsp_);
    s.pack("Blocksqp::exact_hess_lag_sp_", exact_hess_lag_sp_);
    s.pack("Blocksqp::linsol_plugin", linsol_plugin_);
    s.pack("Blocksqp::print_header", print_header_);
    s.pack("Blocksqp::print_iteration", print_iteration_);
    s.pack("Blocksqp::eps", eps_);
    s.pack("Blocksqp::opttol", opttol_);
    s.pack("Blocksqp::nlinfeastol", nlinfeastol_);
    s.pack("Blocksqp::schur", schur_);
    s.pack("Blocksqp::globalization", globalization_);
    s.pack("Blocksqp::restore_feas", restore_feas_);
    s.pack("Blocksqp::max_line_search", max_line_search_);
    s.pack("Blocksqp::max_consec_reduced_steps", max_consec_reduced_steps_);
    s.pack("Blocksqp::max_consec_skipped_updates", max_consec_skipped_updates_);
    s.pack("Blocksqp::max_it_qp", max_it_qp_);
    s.pack("Blocksqp::max_iter", max_iter_);
    s.pack("Blocksqp::warmstart", warmstart_);
    s.pack("Blocksqp::qp_init", qp_init_);
    s.pack("Blocksqp::block_hess", block_hess_);
    s.pack("Blocksqp::hess_scaling", hess_scaling_);
    s.pack("Blocksqp::fallback_scaling", fallback_scaling_);
    s.pack("Blocksqp::max_time_qp", max_time_qp_);
    s.pack("Blocksqp::ini_hess_diag", ini_hess_diag_);
    s.pack("Blocksqp::col_eps", col_eps_);
    s.pack("Blocksqp::col_tau1", col_tau1_);
    s.pack("Blocksqp::col_tau2", col_tau2_);
    s.pack("Blocksqp::hess_damp", hess_damp_);
    s.pack("Blocksqp::hess_damp_fac", hess_damp_fac_);
    s.pack("Blocksqp::hess_update", hess_update_);
    s.pack("Blocksqp::fallback_update", fallback_update_);
    s.pack("Blocksqp::hess_lim_mem", hess_lim_mem_);
    s.pack("Blocksqp::hess_memsize", hess_memsize_);
    s.pack("Blocksqp::which_second_derv", which_second_derv_);
    s.pack("Blocksqp::skip_first_globalization", skip_first_globalization_);
    s.pack("Blocksqp::conv_strategy", conv_strategy_);
    s.pack("Blocksqp::max_conv_qp", max_conv_qp_);
    s.pack("Blocksqp::max_soc_iter", max_soc_iter_);
    s.pack("Blocksqp::gamma_theta", gamma_theta_);
    s.pack("Blocksqp::gamma_f", gamma_f_);
    s.pack("Blocksqp::kappa_soc", kappa_soc_);
    s.pack("Blocksqp::kappa_f", kappa_f_);
    s.pack("Blocksqp::theta_max", theta_max_);
    s.pack("Blocksqp::theta_min", theta_min_);
    s.pack("Blocksqp::delta", delta_);
    s.pack("Blocksqp::s_theta", s_theta_);
    s.pack("Blocksqp::s_f", s_f_);
    s.pack("Blocksqp::kappa_minus", kappa_minus_);
    s.pack("Blocksqp::kappa_plus", kappa_plus_);
    s.pack("Blocksqp::kappa_plus_max", kappa_plus_max_);
    s.pack("Blocksqp::delta_h0", delta_h0_);
    s.pack("Blocksqp::eta", eta_);
    s.pack("Blocksqp::obj_lo", obj_lo_);
    s.pack("Blocksqp::obj_up", obj_up_);
    s.pack("Blocksqp::rho", rho_);
    s.pack("Blocksqp::zeta", zeta_);
    s.pack("Blocksqp::rp_solver", rp_solver_);
    s.pack("Blocksqp::print_maxit_reached", print_maxit_reached_);
  }

} // namespace casadi
