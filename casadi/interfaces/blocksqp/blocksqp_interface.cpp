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


#include "blocksqp_interface.hpp"
#include "casadi/core/std_vector_tools.hpp"

using namespace std;
namespace casadi {

  BlocksqpProblem::BlocksqpProblem(const Blocksqp& self, BlocksqpMemory* m)
    : self(self), m(m) {
    nVar = self.nx_;
    nCon = self.ng_;

    // Bounds on variables and constraints
    bl.Dimension(nVar + nCon).Initialize(-inf);
    bu.Dimension(nVar + nCon).Initialize(inf);
    for (int i=0; i<self.nx_; ++i) {
      bl(i) = m->lbx ? m->lbx[i] : 0;
      bu(i) = m->ubx ? m->ubx[i] : 0;
    }
    for (int i=0; i<self.ng_; ++i) {
      bl(self.nx_ + i) = m->lbg ? m->lbg[i] : 0;
      bu(self.nx_ + i) = m->ubg ? m->ubg[i] : 0;
    }

    // Bounds on objective function
    objLo = -inf;
    objUp = inf;
  }

  void BlocksqpProblem::initialize(blocksqp::Matrix &xi, blocksqp::Matrix &lambda,
                                   blocksqp::Matrix &constrJac) {
    casadi_error("BlocksqpProblem::initialize (dense)");
  }

  void BlocksqpProblem::initialize(blocksqp::Matrix &xi, blocksqp::Matrix &lambda,
                                   double *&jacNz, int *&jacIndRow, int *&jacIndCol) {
    // Primal-dual initial guess
    double* x = xi.array;
    double* lam_x = lambda.array;
    double* lam_g = lam_x + self.nx_;
    casadi_copy(m->x0, self.nx_, x);
    casadi_copy(m->lam_x0, self.nx_, lam_x);
    casadi_copy(m->lam_g0, self.ng_, lam_g);

    // Get Jacobian sparsity
    jacIndRow = new int[self.sp_jac_.nnz()];
    copy_n(self.sp_jac_.row(), self.sp_jac_.nnz(), jacIndRow);
    jacIndCol = const_cast<int*>(self.sp_jac_.colind());
    jacNz = new double[self.sp_jac_.nnz()];
  }

  void BlocksqpProblem::evaluate(const blocksqp::Matrix &xi, const blocksqp::Matrix &lambda,
                                 double *objval, blocksqp::Matrix &constr,
                                 blocksqp::Matrix &gradObj, double *&jacNz, int *&jacIndRow,
                                 int *&jacIndCol,
                                 blocksqp::SymMatrix *&hess, int dmode, int *info) {
    if (dmode==0) {
      // No derivatives
      m->arg[0] = xi.array; // x
      m->arg[1] = m->p; // p
      m->res[0] = objval; // f
      m->res[1] = constr.array; // g
      self.calc_function(m, "nlp_fg");
    } else if (dmode==1) {
      // First order derivatives
      m->arg[0] = xi.array; // x
      m->arg[1] = m->p; // p
      m->res[0] = objval; // f
      m->res[1] = constr.array; // g
      m->res[2] = gradObj.array; // grad:f:x
      m->res[3] = jacNz; // jac:g:x
      self.calc_function(m, "nlp_gf_jg");
    } else {
      casadi_error("Not implemented");
    }
    *info = 0;
  }

  void BlocksqpProblem::evaluate(const blocksqp::Matrix &xi, const blocksqp::Matrix &lambda,
                                 double *objval, blocksqp::Matrix &constr,
                                 blocksqp::Matrix &gradObj, blocksqp::Matrix &constrJac,
                                 blocksqp::SymMatrix *&hess,
                                 int dmode, int *info) {
    if (dmode==0) {
      double *jacNz = 0;
      int *jacIndRow = 0;
      int *jacIndCol = 0;
      return evaluate(xi, lambda, objval, constr, gradObj, jacNz, jacIndRow,
                      jacIndCol, hess, dmode, info);
    }

    casadi_error("BlocksqpProblem::evaluate (dense)");
  }

  void BlocksqpProblem::evaluate(const blocksqp::Matrix &xi, double *objval,
                                 blocksqp::Matrix &constr, int *info) {
    blocksqp::Matrix lambdaDummy, gradObjDummy;
    blocksqp::SymMatrix *hessDummy;
    int dmode = 0;

    blocksqp::Matrix constrJacDummy;
    double *jacNzDummy;
    int *jacIndRowDummy, *jacIndColDummy;
    *info = 0;

    // Try sparse version first
    evaluate(xi, lambdaDummy, objval, constr, gradObjDummy, jacNzDummy,
      jacIndRowDummy, jacIndColDummy, hessDummy, dmode, info);

    // If sparse version is not implemented, try dense version
    if (info) {
      evaluate(xi, lambdaDummy, objval, constr, gradObjDummy,
        constrJacDummy, hessDummy, dmode, info);
    }
  }

  extern "C"
  int CASADI_NLPSOL_BLOCKSQP_EXPORT
  casadi_register_nlpsol_blocksqp(Nlpsol::Plugin* plugin) {
    plugin->creator = Blocksqp::creator;
    plugin->name = "blocksqp";
    plugin->doc = Blocksqp::meta_doc.c_str();
    plugin->version = 30;
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
    clear_memory();
  }

  Options Blocksqp::options_
  = {{&Nlpsol::options_},
     {{"print_level",
       {OT_INT,
        "Print level"}},
      {"debug_level",
       {OT_INT,
        "Amount of debug information that is printed during every iteration"}},
      {"eps",
       {OT_DOUBLE,
        "Values smaller than this are regarded as numerically zero"}},
      {"opttol",
       {OT_DOUBLE,
        "Optimality tolerance"}},
      {"nlinfeastol",
       {OT_DOUBLE,
        "Nonlinear feasibility tolerance"}},
      {"sparse_qp",
       {OT_INT,
        "Which qpOASES variant is used (dense/sparse/Schur)"}},
      {"globalization",
       {OT_INT,
        "Globalization strategy"}},
      {"restore_feas",
       {OT_INT,
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
        "Filter line search parameter, cf. IPOPT paper"}}
     }
  };

  void Blocksqp::init(const Dict& opts) {
    // Call the init method of the base class
    Nlpsol::init(opts);

    // Set default options
    print_level_ = 2;
    debug_level_ = 0;
    eps_ = 1.0e-16;
    opttol_ = 1.0e-6;
    nlinfeastol_ = 1.0e-6;
    sparse_qp_ = 2;
    globalization_ = 1;
    restore_feas_ = 1;
    max_line_search_ = 20;
    max_consec_reduced_steps_ = 100;
    max_consec_skipped_updates_ = 100;
    max_it_qp_ = 5000;
    block_hess_ = 1;
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

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="print_level") {
        print_level_ = op.second;
      } else if (op.first=="debug_level") {
        debug_level_ = op.second;
      } else if (op.first=="eps") {
        eps_ = op.second;
      } else if (op.first=="opttol") {
        opttol_ = op.second;
      } else if (op.first=="nlinfeastol") {
        nlinfeastol_ = op.second;
      } else if (op.first=="sparse_qp") {
        sparse_qp_ = op.second;
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
      }
    }

    // If we compute second constraints derivatives switch to
    // finite differences Hessian (convenience)
    if (which_second_derv_ == 2) {
      hess_update_ = 4;
      block_hess_ = 1;
    }

    // If we don't use limited memory BFGS we need to store only one vector.
    if (!hess_lim_mem_) hess_memsize_ = 1;
    if (sparse_qp_ != 2 && hess_update_ == 1) {
      printf("SR1 update only works with qpOASES Schur complement version. "
             "Using BFGS updates instead.\n");
      hess_update_ = 2;
      hess_scaling_ = fallback_scaling_;
    }

    // Setup NLP functions
    create_function("nlp_fg", {"x", "p"}, {"f", "g"});
    Function gf_jg = create_function("nlp_gf_jg", {"x", "p"},
                                     {"f", "g", "grad:f:x", "jac:g:x"});
    sp_jac_ = gf_jg.sparsity_out("jac_g_x");

    if (block_hess_ == 0) {
      // No block-structured Hessian
      blocks_ = {0, nx_};
      which_second_derv_ = 0;
    } else {
      // Detect block structure

      // Get the sparsity pattern for the Hessian of the Lagrangian
      Function grad_lag = oracle_.factory("grad_lag",
                                          {"x", "p", "lam:f", "lam:g"}, {"grad:gamma:x"},
                                          {{"gamma", {"f", "g"}}});
      Sparsity Hsp = grad_lag.sparsity_jac("x", "grad_gamma_x", false, true);

      // Make sure diagonal exists
      Hsp = Hsp + Sparsity::diag(nx_);

      // Find the strongly connected components of the Hessian
      // Unlike Sparsity::scc, assume ordered
      const int* colind = Hsp.colind();
      const int* row = Hsp.row();
      blocks_.push_back(0);
      int ind = 0;
      while (ind < nx_) {
        // Find the next cutoff
        int next=ind+1;
        while (ind<next && ind<nx_) {
          for (int k=colind[ind]; k<colind[ind+1]; ++k) next = max(next, 1+row[k]);
          ind++;
        }
        blocks_.push_back(next);
      }

      // hybrid strategy: 1 block for constraints, 1 for objective
      if (block_hess_ == 2 && blocks_.size() > 3) {
        blocks_ = {0, *(blocks_.rbegin()+1), blocks_.back()};
      }
    }

    // Number of blocks
    nblocks_ = blocks_.size()-1;

    // Allocate memory
    alloc_w(sp_jac_.nnz(), true); // jac
  }

  void Blocksqp::init_memory(void* mem) const {
    Nlpsol::init_memory(mem);
    auto m = static_cast<BlocksqpMemory*>(mem);
  }

  void Blocksqp::set_work(void* mem, const double**& arg, double**& res,
                                   int*& iw, double*& w) const {
    auto m = static_cast<BlocksqpMemory*>(mem);

    // Set work in base classes
    Nlpsol::set_work(mem, arg, res, iw, w);

    // Temporary memory
    m->jac = w; w += sp_jac_.nnz();
  }

  void Blocksqp::solve(void* mem) const {
    auto m = static_cast<BlocksqpMemory*>(mem);

    int ret = 0;
    char outpath[255];
    strcpy(outpath, "./");

    // Create problem evaluation object
    vector<int> blocks = blocks_;
    m->prob = new BlocksqpProblem(*this, m);

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

    int maxblocksize = 1;

    for (int k=0; k<nblocks_+1; k++) {
      if (k > 0)
        if (blocks_[k] - blocks_[k-1] > maxblocksize)
          maxblocksize = blocks_[k] - blocks_[k-1];
    }

    if (hess_lim_mem_ && hess_memsize_ == 0)
      const_cast<Blocksqp*>(this)->hess_memsize_ = maxblocksize;

    allocMin(m, m->prob);

    if (!sparse_qp_) {
      m->constrJac.Dimension(m->prob->nCon, m->prob->nVar).Initialize(0.0);
      m->hessNz = new double[m->prob->nVar * m->prob->nVar];
    } else {
      m->hessNz = 0;
    }

    m->jacNz = 0;
    m->jacIndCol = 0;
    m->jacIndRow = 0;

    m->hessIndCol = 0;
    m->hessIndRow = 0;
    m->hessIndLo = 0;
    m->hess = 0;
    m->hess1 = 0;
    m->hess2 = 0;

    m->noUpdateCounter = 0;

    allocHess(m);
    allocAlg(m, m->prob);

    if (sparse_qp_ < 2) {
      m->qp = new qpOASES::SQProblem(m->prob->nVar, m->prob->nCon);
      m->qpSave = new qpOASES::SQProblem(m->prob->nVar, m->prob->nCon);
    } else {
      m->qp = new qpOASES::SQProblemSchur(m->prob->nVar, m->prob->nCon, qpOASES::HST_UNKNOWN, 50);
      m->qpSave = new qpOASES::SQProblemSchur(m->prob->nVar,
        m->prob->nCon, qpOASES::HST_UNKNOWN, 50);
    }

    m->initCalled = false;

    // Print header and information about the algorithmic parameters
    printInfo(m, print_level_);

    // Open output files
    initStats(m);
    initIterate(m);

    // Initialize filter with pair (maxConstrViolation, objLowerBound)
    initializeFilter(m);

    // Set initial values for all xi and set the Jacobian for linear constraints
    if (sparse_qp_) {
      m->prob->initialize(m->xi, m->lambda, m->jacNz,
        m->jacIndRow, m->jacIndCol);
    } else {
      m->prob->initialize(m->xi, m->lambda, m->constrJac);
    }

    m->initCalled = true;




    ret = run(m, 100);
    finish(m);
    if (ret==1) casadi_warning("Maximum number of iterations reached");

    // Get optimal cost
    if (m->f) *m->f = m->obj;
    // Get primal solution
    casadi_copy(m->xi.array, nx_, m->x);
    // Get dual solution (simple bounds)
    if (m->lam_x) {
      casadi_copy(m->lambda.array, nx_, m->lam_x);
      casadi_scal(nx_, -1., m->lam_x);
    }
    // Get dual solution (nonlinear bounds)
    if (m->lam_g) {
      casadi_copy(m->lambda.array + nx_, ng_, m->lam_g);
      casadi_scal(ng_, -1., m->lam_g);
    }

    // Clean up
    delete m->prob;
    delete m->qp;
    delete m->qpSave;
    if (m->noUpdateCounter != 0) delete[] m->noUpdateCounter;
    if (m->jacNz != 0) delete[] m->jacNz;
    if (m->jacIndRow != 0) delete[] m->jacIndRow;
    if (m->hessNz != 0) delete[] m->hessNz;
    if (m->hessIndRow != 0) delete[] m->hessIndRow;
  }

  int Blocksqp::run(BlocksqpMemory* m, int maxIt, int warmStart) const {
    int it, infoQP = 0, infoEval = 0;
    bool skipLineSearch = false;
    bool hasConverged = false;
    int whichDerv = which_second_derv_;

    if (!m->initCalled) {
      printf("init() must be called before run(). Aborting.\n");
      return -1;
    }

    if (warmStart == 0 || m->itCount == 0) {
      // SQP iteration 0

      /// Set initial Hessian approximation
      calcInitialHessian(m);

      /// Evaluate all functions and gradients for xi_0
      if (sparse_qp_) {
        m->prob->evaluate(m->xi, m->lambda, &m->obj,
          m->constr, m->gradObj,
                        m->jacNz, m->jacIndRow, m->jacIndCol,
                        m->hess, 1+whichDerv, &infoEval);
      } else {
        m->prob->evaluate(m->xi, m->lambda, &m->obj,
          m->constr, m->gradObj,
                        m->constrJac, m->hess, 1+whichDerv, &infoEval);
      }
      m->nDerCalls++;

      /// Check if converged
      hasConverged = calcOptTol(m);
      printProgress(m, m->prob, hasConverged);
      if (hasConverged)
        return 0;

      m->itCount++;
    }

    /*
     * SQP Loop: during first iteration, m->itCount = 1
     */
    for (it=0; it<maxIt; it++) {
      /// Solve QP subproblem with qpOASES or QPOPT
      updateStepBounds(m, 0);
      infoQP = solveQP(m, m->deltaXi, m->lambdaQP);

      if (infoQP == 1) {
          // 1.) Maximum number of iterations reached
          printf("***Warning! Maximum number of QP iterations exceeded.***\n");
      } else if (infoQP == 2 || infoQP > 3) {
          // 2.) QP error (e.g., unbounded), solve again with pos.def. diagonal matrix (identity)
          printf("***QP error. Solve again with identity matrix.***\n");
          resetHessian(m);
          infoQP = solveQP(m, m->deltaXi, m->lambdaQP);
          if (infoQP) {
            // If there is still an error, terminate.
            printf("***QP error. Stop.***\n");
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
          printf("***QP infeasible. Trying to reduce constraint violation...");
          qpError = feasibilityRestorationHeuristic(m);
          if (!qpError) {
            m->steptype = 2;
            printf("Success.***\n");
          } else {
            printf("Failed.***\n");
          }
        }

        // Invoke feasibility restoration phase
        //if (qpError && m->steptype < 3 && restore_feas_)
        if (qpError && restore_feas_ && m->cNorm > 0.01 * nlinfeastol_) {
          printf("***Start feasibility restoration phase.***\n");
          m->steptype = 3;
          qpError = feasibilityRestorationPhase(m);
        }

        // If everything failed, abort.
        if (qpError) {
          printf("***QP error. Stop.***\n");
          return -1;
        }
      }

      /// Determine steplength alpha
      if (globalization_ == 0 || (skip_first_globalization_
        && m->itCount == 1)) {
        // No globalization strategy, but reduce step if function cannot be evaluated
        if (fullstep(m)) {
          printf("***Constraint or objective could not be evaluated at new point. Stop.***\n");
          return -1;
        }
        m->steptype = 0;
      } else if (globalization_ == 1 && !skipLineSearch) {
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

            printf("***Warning! Steplength too short. Trying to reduce constraint violation...");

            // Integration over whole time interval
            lsError = feasibilityRestorationHeuristic(m);
            if (!lsError) {
                m->steptype = 2;
                printf("Success.***\n");
              } else {
              printf("Failed.***\n");
            }
          }

          // Heuristic 3: Recompute step with a diagonal Hessian
          if (lsError && m->steptype != 1 && m->steptype != 2) {
            // After closing continuity gaps, we already take a step with initial Hessian.
            // If this step is not accepted then this will cause an infinite loop!

            printf("***Warning! Steplength too short. "
                  "Trying to find a new step with identity Hessian.***\n");
            m->steptype = 1;

            resetHessian(m);
            continue;
          }

          // If this does not yield a successful step, start restoration phase
          if (lsError && m->cNorm > 0.01 * nlinfeastol_ && restore_feas_) {
            printf("***Warning! Steplength too short. Start feasibility restoration phase.***\n");
            m->steptype = 3;

            // Solve NLP with minimum norm objective
            lsError = feasibilityRestorationPhase(m);
          }

          // If everything failed, abort.
          if (lsError) {
            printf("***Line search error. Stop.***\n");
            return -1;
          }
        } else {
          m->steptype = 0;
        }
      }

      /// Calculate "old" Lagrange gradient: gamma = dL(xi_k, lambda_k+1)
      calcLagrangeGradient(m, m->gamma, 0);

      /// Evaluate functions and gradients at the new xi
      if (sparse_qp_) {
        m->prob->evaluate(m->xi, m->lambda, &m->obj, m->constr,
          m->gradObj, m->jacNz, m->jacIndRow,
          m->jacIndCol, m->hess, 1+whichDerv, &infoEval);
      } else {
        m->prob->evaluate(m->xi, m->lambda, &m->obj, m->constr,
            m->gradObj, m->constrJac, m->hess, 1+whichDerv, &infoEval);
      }
      m->nDerCalls++;

      /// Check if converged
      hasConverged = calcOptTol(m);

      /// Print one line of output for the current iteration
      printProgress(m, m->prob, hasConverged);
      if (hasConverged && m->steptype < 2) {
        m->itCount++;
        if (debug_level_ > 2) {
          //printf("Computing finite differences Hessian at the solution ... \n");
          //calcFiniteDiffHessian();
          //m->printHessian(nblocks_, m->hess);
          dumpQPCpp(m, m->prob, m->qp, sparse_qp_);
        }
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
        casadi_error("Not implemented");
      }

      // If limited memory updates  are used, set pointers deltaXi and
      // gamma to the next column in deltaMat and gammaMat
      updateDeltaGamma(m);

      m->itCount++;
      skipLineSearch = false;
    }

    return 1;
  }


  void Blocksqp::finish(BlocksqpMemory* m) const {
    if (m->initCalled) {
      m->initCalled = false;
    } else {
      printf("init() must be called before finish().\n");
      return;
    }

    if (debug_level_ > 0) {
      fprintf(m->progressFile, "\n");
      fclose(m->progressFile);
      fprintf(m->updateFile, "\n");
      fclose(m->updateFile);
    }

    if (debug_level_ > 1) {
      fclose(m->primalVarsFile);
      fclose(m->dualVarsFile);
    }
  }


  /**
   * Compute gradient of Lagrangian or difference of Lagrangian gradients (sparse version)
   *
   * flag == 0: output dL(xi, lambda)
   * flag == 1: output dL(xi_k+1, lambda_k+1) - L(xi_k, lambda_k+1)
   * flag == 2: output dL(xi_k+1, lambda_k+1) - df(xi)
   */
  void Blocksqp::
  calcLagrangeGradient(BlocksqpMemory* m, const blocksqp::Matrix &lambda,
    const blocksqp::Matrix &gradObj, double *jacNz, int *jacIndRow,
    int *jacIndCol, blocksqp::Matrix &gradLagrange, int flag) const {
    int iVar, iCon;

    // Objective gradient
    if (flag == 0) {
      for (iVar=0; iVar<m->prob->nVar; iVar++) {
        gradLagrange(iVar) = gradObj(iVar);
      }
    } else if (flag == 1) {
      for (iVar=0; iVar<m->prob->nVar; iVar++) {
        gradLagrange(iVar) = gradObj(iVar) - gradLagrange(iVar);
      }
    } else {
      gradLagrange.Initialize(0.0);
    }

    // - lambdaT * constrJac
    for (iVar=0; iVar<m->prob->nVar; iVar++)
      for (iCon=jacIndCol[iVar]; iCon<jacIndCol[iVar+1]; iCon++)
        gradLagrange(iVar) -= lambda(m->prob->nVar + jacIndRow[iCon]) * jacNz[iCon];

    // - lambdaT * simpleBounds
    for (iVar=0; iVar<m->prob->nVar; iVar++) gradLagrange(iVar) -= lambda(iVar);
  }


  /**
   * Compute gradient of Lagrangian or difference of Lagrangian gradients (dense version)
   *
   * flag == 0: output dL(xi, lambda)
   * flag == 1: output dL(xi_k+1, lambda_k+1) - L(xi_k, lambda_k+1)
   * flag == 2: output dL(xi_k+1, lambda_k+1) - df(xi)
   */
  void Blocksqp::
  calcLagrangeGradient(BlocksqpMemory* m, const blocksqp::Matrix &lambda,
                       const blocksqp::Matrix &gradObj, const blocksqp::Matrix &constrJac,
                       blocksqp::Matrix &gradLagrange, int flag) const {
    int iVar, iCon;

    // Objective gradient
    if (flag == 0) {
      for (iVar=0; iVar<m->prob->nVar; iVar++) {
        gradLagrange(iVar) = gradObj(iVar);
      }
    } else if (flag == 1) {
      for (iVar=0; iVar<m->prob->nVar; iVar++) {
        gradLagrange(iVar) = gradObj(iVar) - gradLagrange(iVar);
      }
    } else {
      gradLagrange.Initialize(0.0);
    }

    // - lambdaT * constrJac
    for (iVar=0; iVar<m->prob->nVar; iVar++)
      for (iCon=0; iCon<m->prob->nCon; iCon++)
        gradLagrange(iVar) -= lambda(m->prob->nVar + iCon) * constrJac(iCon, iVar);

    // - lambdaT * simpleBounds
    for (iVar=0; iVar<m->prob->nVar; iVar++) {
      gradLagrange(iVar) -= lambda(iVar);
    }
  }


  /**
   * Wrapper if called with standard arguments
   */
  void Blocksqp::
  calcLagrangeGradient(BlocksqpMemory* m, blocksqp::Matrix &gradLagrange, int flag) const {
    if (sparse_qp_) {
      calcLagrangeGradient(m, m->lambda, m->gradObj, m->jacNz,
        m->jacIndRow, m->jacIndCol, gradLagrange, flag);
    } else {
      calcLagrangeGradient(m, m->lambda, m->gradObj, m->constrJac,
        gradLagrange, flag);
    }
  }


  /**
   * Compute optimality conditions:
   * ||gradLagrange(xi,lambda)||_infty / (1 + ||lambda||_infty) <= TOL
   * and
   * ||constrViolation||_infty / (1 + ||xi||_infty) <= TOL
   */
  bool Blocksqp::calcOptTol(BlocksqpMemory* m) const {
    // scaled norm of Lagrangian gradient
    calcLagrangeGradient(m, m->gradLagrange, 0);
    m->gradNorm = lInfVectorNorm(m->gradLagrange);
    m->tol = m->gradNorm /(1.0 + lInfVectorNorm(m->lambda));

    // norm of constraint violation
    m->cNorm  = lInfConstraintNorm(m->xi, m->constr, m->prob->bu, m->prob->bl);
    m->cNormS = m->cNorm /(1.0 + lInfVectorNorm(m->xi));

    if (m->tol <= opttol_ && m->cNormS <= nlinfeastol_)
      return true;
    else
      return false;
  }

  void Blocksqp::printInfo(BlocksqpMemory* m, int printLevel) const {
    char hessString1[100];
    char hessString2[100];
    char globString[100];
    char qpString[100];

    if (printLevel == 0)
      return;

    /* QP Solver */
    if (sparse_qp_ == 0)
      strcpy(qpString, "dense, reduced Hessian factorization");
    else if (sparse_qp_ == 1)
      strcpy(qpString, "sparse, reduced Hessian factorization");
    else if (sparse_qp_ == 2)
      strcpy(qpString, "sparse, Schur complement approach");

    /* Globalization */
    if (globalization_ == 0)
      strcpy(globString, "none (full step)");
    else if (globalization_ == 1)
      strcpy(globString, "filter line search");

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

    printf("\n+---------------------------------------------------------------+\n");
    printf("| Starting blockSQP with the following algorithmic settings:    |\n");
    printf("+---------------------------------------------------------------+\n");
    printf("| qpOASES flavor            | %-34s|\n", qpString);
    printf("| Globalization             | %-34s|\n", globString);
    printf("| 1st Hessian approximation | %-34s|\n", hessString1);
    printf("| 2nd Hessian approximation | %-34s|\n", hessString2);
    printf("+---------------------------------------------------------------+\n\n");
  }

  void Blocksqp::
  acceptStep(BlocksqpMemory* m, double alpha) const {
    acceptStep(m, m->deltaXi, m->lambdaQP, alpha, 0);
  }

  void Blocksqp::
  acceptStep(BlocksqpMemory* m, const blocksqp::Matrix &deltaXi,
    const blocksqp::Matrix &lambdaQP, double alpha, int nSOCS) const {
    int k;
    double lStpNorm;

    // Current alpha
    m->alpha = alpha;
    m->nSOCS = nSOCS;

    // Set new xi by accepting the current trial step
    for (k=0; k<m->xi.M(); k++) {
      m->xi(k) = m->trialXi(k);
      m->deltaXi(k) = alpha * deltaXi(k);
    }

    // Store the infinity norm of the multiplier step
    m->lambdaStepNorm = 0.0;
    for (k=0; k<m->lambda.M(); k++)
      if ((lStpNorm = fabs(alpha*lambdaQP(k) - alpha*m->lambda(k))) > m->lambdaStepNorm)
        m->lambdaStepNorm = lStpNorm;

    // Set new multipliers
    for (k=0; k<m->lambda.M(); k++)
      m->lambda(k) = (1.0 - alpha)*m->lambda(k) + alpha*lambdaQP(k);

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
    int i;
    int nVar = m->prob->nVar;

    // Update bounds on linearized constraints for the next SOC QP:
    // That is different from the update for the first SOC QP!
    for (i=0; i<m->prob->nCon; i++) {
      if (m->prob->bl(nVar+i) != inf)
        m->deltaBl(nVar+i) = (*alphaSOC)*m->deltaBl(nVar+i) - m->constr(i);
      else
        m->deltaBl(nVar+i) = inf;

      if (m->prob->bu(nVar+i) != inf)
        m->deltaBu(nVar+i) = (*alphaSOC)*m->deltaBu(nVar+i) - m->constr(i);
      else
        m->deltaBu(nVar+i) = inf;
    }

    *alphaSOC = (*alphaSOC) * 0.5;
  }


  /**
   * Take a full Quasi-Newton step, except when integrator fails:
   * xi = xi + deltaXi
   * lambda = lambdaQP
   */
  int Blocksqp::fullstep(BlocksqpMemory* m) const {
    double alpha;
    double objTrial, cNormTrial;
    int i, k, info;
    int nVar = m->prob->nVar;

    // Backtracking line search
    alpha = 1.0;
    for (k=0; k<10; k++) {
      // Compute new trial point
      for (i=0; i<nVar; i++)
        m->trialXi(i) = m->xi(i) + alpha * m->deltaXi(i);

      // Compute problem functions at trial point
      m->prob->evaluate(m->trialXi, &objTrial, m->constr, &info);
      m->nFunCalls++;
      cNormTrial = lInfConstraintNorm(m->trialXi, m->constr, m->prob->bu, m->prob->bl);
      // Reduce step if evaluation fails, if lower bound is violated
      // or if objective or a constraint is NaN
      if (info != 0 || objTrial < m->prob->objLo || objTrial > m->prob->objUp
        || !(objTrial == objTrial) || !(cNormTrial == cNormTrial)) {
        printf("info=%i, objTrial=%g\n", info, objTrial);
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
  int Blocksqp::filterLineSearch(BlocksqpMemory* m) const {
    double alpha = 1.0;
    double cNorm, cNormTrial, objTrial, dfTdeltaXi;

    int i, k, info;
    int nVar = m->prob->nVar;

    // Compute ||constr(xi)|| at old point
    cNorm = lInfConstraintNorm(m->xi, m->constr, m->prob->bu, m->prob->bl);

    // Backtracking line search
    for (k=0; k<max_line_search_; k++) {
      // Compute new trial point
      for (i=0; i<nVar; i++)
        m->trialXi(i) = m->xi(i) + alpha * m->deltaXi(i);

      // Compute grad(f)^T * deltaXi
      dfTdeltaXi = 0.0;
      for (i=0; i<nVar; i++)
        dfTdeltaXi += m->gradObj(i) * m->deltaXi(i);

      // Compute objective and at ||constr(trialXi)||_1 at trial point
      m->prob->evaluate(m->trialXi, &objTrial, m->constr, &info);
      m->nFunCalls++;
      cNormTrial = lInfConstraintNorm(m->trialXi, m->constr, m->prob->bu, m->prob->bl);
      // Reduce step if evaluation fails, if lower bound is violated or if objective is NaN
      if (info != 0 || objTrial < m->prob->objLo || objTrial > m->prob->objUp
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
   * min_d d^TBd + d^TgradObj
   * s.t.  bl <= A^Td + constr(xi+alpha*deltaXi) - A^TdeltaXi <= bu
   *
   */
  bool Blocksqp::
  secondOrderCorrection(BlocksqpMemory* m, double cNorm, double cNormTrial,
    double dfTdeltaXi, bool swCond, int it) const {

    // Perform SOC only on the first iteration of backtracking line search
    if (it > 0) return false;
    // If constraint violation of the trialstep is lower than the current one skip SOC
    if (cNormTrial < cNorm) return false;

    int nSOCS = 0;
    double cNormTrialSOC, cNormOld, objTrialSOC;
    int i, k, info;
    int nVar = m->prob->nVar;
    blocksqp::Matrix deltaXiSOC, lambdaQPSOC;

    // m->constr contains result at first trial point: c(xi+deltaXi)
    // m->constrJac, m->AdeltaXi and m->gradObj are unchanged so far.

    // First SOC step
    deltaXiSOC.Dimension(m->deltaXi.M()).Initialize(0.0);
    lambdaQPSOC.Dimension(m->lambdaQP.M()).Initialize(0.0);

    // Second order correction loop
    cNormOld = cNorm;
    for (k=0; k<max_soc_iter_; k++) {
      nSOCS++;

      // Update bounds for SOC QP
      updateStepBounds(m, 1);

      // Solve SOC QP to obtain new, corrected deltaXi
      // (store in separate vector to avoid conflict with original deltaXi
      // -> need it in linesearch!)
      info = solveQP(m, deltaXiSOC, lambdaQPSOC, false);
      if (info != 0) return false; // Could not solve QP, abort SOC

      // Set new SOC trial point
      for (i=0; i<nVar; i++) {
        m->trialXi(i) = m->xi(i) + deltaXiSOC(i);
      }

      // Compute objective and ||constr(trialXiSOC)||_1 at SOC trial point
      m->prob->evaluate(m->trialXi, &objTrialSOC, m->constr, &info);
      m->nFunCalls++;
      cNormTrialSOC = lInfConstraintNorm(m->trialXi, m->constr,
        m->prob->bu, m->prob->bl);
      if (info != 0 || objTrialSOC < m->prob->objLo || objTrialSOC > m->prob->objUp
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
          acceptStep(m, deltaXiSOC, lambdaQPSOC, 1.0, nSOCS);
          return true;
        }
      }

      // Check sufficient decrease, case II (in SOC)
      if (cNorm > theta_min_ || !swCond) {
        if (cNormTrialSOC < (1.0 - gamma_theta_) * cNorm
        || objTrialSOC < m->obj - gamma_f_ * cNorm) {
          // found suitable alpha during SOC, stop
          acceptStep(m, deltaXiSOC, lambdaQPSOC, 1.0, nSOCS);
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
  int Blocksqp::feasibilityRestorationPhase(BlocksqpMemory* m) const {
    // No Feasibility restoration phase
    if (restore_feas_ == 0) return -1;

    casadi_error("not implemented");
    return 0;
  }


  /**
   * Try to (partly) improve constraint violation by satisfying
   * the (pseudo) continuity constraints, i.e. do a single shooting
   * iteration with the current controls and measurement weights q and w
   */
  int Blocksqp::feasibilityRestorationHeuristic(BlocksqpMemory* m) const {
    m->nRestHeurCalls++;

    int info, k;
    double cNormTrial;

    info = 0;

    // Call problem specific heuristic to reduce constraint violation.
    // For shooting methods that means setting consistent values for
    // shooting nodes by one forward integration.
    for (k=0; k<m->prob->nVar; k++) // input: last successful step
      m->trialXi(k) = m->xi(k);
    m->prob->reduceConstrVio(m->trialXi, &info);
    if (info) {
      // If an error occured in restoration heuristics, abort
      return -1;
    }

    // Compute objective and constraints at the new (hopefully feasible) point
    m->prob->evaluate(m->trialXi, &m->obj, m->constr, &info);
    m->nFunCalls++;
    cNormTrial = lInfConstraintNorm(m->trialXi, m->constr, m->prob->bu, m->prob->bl);
    if (info != 0 || m->obj < m->prob->objLo || m->obj > m->prob->objUp
      || !(m->obj == m->obj) || !(cNormTrial == cNormTrial))
      return -1;

    // Is the new point acceptable for the filter?
    if (pairInFilter(m, cNormTrial, m->obj)) {
      // point is in the taboo region, restoration heuristic not successful!
      return -1;
    }

    // If no error occured in the integration all shooting variables now
    // have the values obtained by a single shooting integration.
    // This is done instead of a Newton-like step in the current SQP iteration

    m->alpha = 1.0;
    m->nSOCS = 0;

    // reset reduced step counter
    m->reducedStepCount = 0;

    // Reset lambda
    m->lambda.Initialize(0.0);
    m->lambdaQP.Initialize(0.0);

    // Compute the "step" taken by closing the continuity conditions
    /// \note deltaXi is reset by resetHessian(), so this doesn't matter
    for (k=0; k<m->prob->nVar; k++) {
      //m->deltaXi(k) = m->trialXi(k) - m->xi(k);
      m->xi(k) = m->trialXi(k);
    }

    // reduce Hessian and limited memory information
    resetHessian(m);

    return 0;
  }


  /**
   * If the line search fails, check if the full step reduces the KKT error by a factor kappaF.
   */
  int Blocksqp::kktErrorReduction(BlocksqpMemory* m) const {
    int i, info = 0;
    double objTrial, cNormTrial, trialGradNorm, trialTol;
    blocksqp::Matrix trialConstr, trialGradLagrange;

    // Compute new trial point
    for (i=0; i<m->prob->nVar; i++)
      m->trialXi(i) = m->xi(i) + m->deltaXi(i);

    // Compute objective and ||constr(trialXi)|| at trial point
    trialConstr.Dimension(m->prob->nCon).Initialize(0.0);
    m->prob->evaluate(m->trialXi, &objTrial, trialConstr, &info);
    m->nFunCalls++;
    cNormTrial = lInfConstraintNorm(m->trialXi, trialConstr, m->prob->bu, m->prob->bl);
    if (info != 0 || objTrial < m->prob->objLo || objTrial > m->prob->objUp
      || !(objTrial == objTrial) || !(cNormTrial == cNormTrial)) {
      // evaluation error
      return 1;
    }

    // Compute KKT error of the new point

    // scaled norm of Lagrangian gradient
    trialGradLagrange.Dimension(m->prob->nVar).Initialize(0.0);
    if (sparse_qp_) {
      calcLagrangeGradient(m, m->lambdaQP, m->gradObj, m->jacNz,
                            m->jacIndRow, m->jacIndCol, trialGradLagrange, 0);
    } else {
      calcLagrangeGradient(m, m->lambdaQP, m->gradObj, m->constrJac,
                            trialGradLagrange, 0);
    }

    trialGradNorm = lInfVectorNorm(trialGradLagrange);
    trialTol = trialGradNorm /(1.0 + lInfVectorNorm(m->lambdaQP));

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
    std::set< std::pair<double, double> >::iterator iter;
    std::set< std::pair<double, double> > *filter;
    filter = m->filter;

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

    for (iter=filter->begin(); iter!=filter->end(); iter++)
      if ((cNorm >= (1.0 - gamma_theta_) * iter->first ||
           (cNorm < 0.01 * nlinfeastol_ && iter->first < 0.01 * nlinfeastol_)) &&
          obj >= iter->second - gamma_f_ * iter->first) {
        return 1;
      }

    return 0;
  }


  void Blocksqp::initializeFilter(BlocksqpMemory* m) const {
    std::set< std::pair<double, double> >::iterator iter;
    std::pair<double, double> initPair(theta_max_, m->prob->objLo);

    // Remove all elements
    iter=m->filter->begin();
    while (iter != m->filter->end()) {
      std::set< std::pair<double, double> >::iterator iterToRemove = iter;
      iter++;
      m->filter->erase(iterToRemove);
    }

    // Initialize with pair (maxConstrViolation, objLowerBound);
    m->filter->insert(initPair);
  }


  /**
   * Augment the filter:
   * F_k+1 = F_k U { (c,f) | c > (1-gammaTheta)cNorm and f > obj-gammaF*c
   */
  void Blocksqp::
  augmentFilter(BlocksqpMemory* m, double cNorm, double obj) const {
    std::set< std::pair<double, double> >::iterator iter;
    std::pair<double, double> entry((1.0 - gamma_theta_)*cNorm, obj
      - gamma_f_*cNorm);

    // Augment filter by current element
    m->filter->insert(entry);

    // Remove dominated elements
    iter=m->filter->begin();
    while (iter != m->filter->end()) {
      if (iter->first > entry.first && iter->second > entry.second) {
        std::set< std::pair<double, double> >::iterator iterToRemove = iter;
        iter++;
        m->filter->erase(iterToRemove);
      } else {
        iter++;
      }
    }
  }

  /**
   * Initial Hessian: Identity matrix
   */
  void Blocksqp::calcInitialHessian(BlocksqpMemory* m) const {
    int iBlock;

    for (iBlock=0; iBlock<nblocks_; iBlock++)
      //if objective derv is computed exactly, don't set the last block!
      if (!(which_second_derv_ == 1 && block_hess_
        && iBlock == nblocks_-1))
        calcInitialHessian(m, iBlock);
  }


  /**
   * Initial Hessian for one block: Identity matrix
   */
  void Blocksqp::calcInitialHessian(BlocksqpMemory* m, int iBlock) const {
    m->hess[iBlock].Initialize(0.0);

    // Each block is a diagonal matrix
    for (int i=0; i<m->hess[iBlock].M(); i++)
      m->hess[iBlock](i, i) = ini_hess_diag_;

    // If we maintain 2 Hessians, also reset the second one
    if (m->hess2 != 0) {
      m->hess2[iBlock].Initialize(0.0);
      for (int i=0; i<m->hess2[iBlock].M(); i++)
        m->hess2[iBlock](i, i) = ini_hess_diag_;
    }
  }


  void Blocksqp::resetHessian(BlocksqpMemory* m) const {
    for (int iBlock=0; iBlock<nblocks_; iBlock++) {
      if (!(which_second_derv_ == 1 && block_hess_ && iBlock == nblocks_ - 1)) {
        // if objective derv is computed exactly, don't set the last block!
        resetHessian(m, iBlock);
      }
    }
  }


  void Blocksqp::resetHessian(BlocksqpMemory* m, int iBlock) const {
    blocksqp::Matrix smallDelta, smallGamma;
    int nVarLocal = m->hess[iBlock].M();

    // smallGamma and smallDelta are either subvectors of gamma and delta
    // or submatrices of gammaMat, deltaMat, i.e. subvectors of gamma and delta
    // from m prev. iterations (for L-BFGS)
    smallGamma.Submatrix(m->gammaMat, nVarLocal, m->gammaMat.N(),
      blocks_[iBlock], 0);
    smallDelta.Submatrix(m->deltaMat, nVarLocal, m->deltaMat.N(),
      blocks_[iBlock], 0);

    // Remove past information on Lagrangian gradient difference
    smallGamma.Initialize(0.0);

    // Remove past information on steps
    smallDelta.Initialize(0.0);

    // Remove information on old scalars (used for COL sizing)
    m->deltaNorm(iBlock) = 1.0;
    m->deltaGamma(iBlock) = 0.0;
    m->deltaNormOld(iBlock) = 1.0;
    m->deltaGammaOld(iBlock) = 0.0;

    m->noUpdateCounter[iBlock] = -1;

    calcInitialHessian(m, iBlock);
  }

  void Blocksqp::sizeInitialHessian(BlocksqpMemory* m, const blocksqp::Matrix &gamma,
    const blocksqp::Matrix &delta, int iBlock, int option) const {
    int i, j;
    double scale;
    double myEps = 1.0e3 * eps_;

    if (option == 1) {
      // Shanno-Phua
      scale = adotb(gamma, gamma) / fmax(adotb(delta, gamma), myEps);
    } else if (option == 2) {
      // Oren-Luenberger
      scale = adotb(delta, gamma) / fmax(adotb(delta, delta), myEps);
      scale = fmin(scale, 1.0);
    } else if (option == 3) {
      // Geometric mean of 1 and 2
      scale = sqrt(adotb(gamma, gamma) / fmax(adotb(delta, delta), myEps));
    } else {
      // Invalid option, ignore
      return;
    }

    if (scale > 0.0) {
      scale = fmax(scale, myEps);
      for (i=0; i<m->hess[iBlock].M(); i++)
        for (j=i; j<m->hess[iBlock].M(); j++)
          m->hess[iBlock](i, j) *= scale;
    } else {
      scale = 1.0;
    }

    // statistics: average sizing factor
    m->averageSizingFactor += scale;
  }


  void Blocksqp::sizeHessianCOL(BlocksqpMemory* m, const blocksqp::Matrix &gamma,
    const blocksqp::Matrix &delta, int iBlock) const {
    int i, j;
    double theta, scale, myEps = 1.0e3 * eps_;
    double deltaNorm, deltaNormOld, deltaGamma, deltaGammaOld, deltaBdelta;

    // Get sTs, sTs_, sTy, sTy_, sTBs
    deltaNorm = m->deltaNorm(iBlock);
    deltaGamma = m->deltaGamma(iBlock);
    deltaNormOld = m->deltaNormOld(iBlock);
    deltaGammaOld = m->deltaGammaOld(iBlock);
    deltaBdelta = 0.0;
    for (i=0; i<delta.M(); i++)
      for (j=0; j<delta.M(); j++)
        deltaBdelta += delta(i) * m->hess[iBlock](i, j) * delta(j);

    // Centered Oren-Luenberger factor
    if (m->noUpdateCounter[iBlock] == -1) {
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
      //printf("Sizing value (COL) block %i = %g\n", iBlock, scale);
      for (i=0; i<m->hess[iBlock].M(); i++)
        for (j=i; j<m->hess[iBlock].M(); j++)
          m->hess[iBlock](i, j) *= scale;

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
  calcHessianUpdate(BlocksqpMemory* m, int updateType, int hessScaling) const {
    int iBlock, nBlocks;
    int nVarLocal;
    blocksqp::Matrix smallGamma, smallDelta;
    bool firstIter;

    //if objective derv is computed exactly, don't set the last block!
    if (which_second_derv_ == 1 && block_hess_)
      nBlocks = nblocks_ - 1;
    else
      nBlocks = nblocks_;

    // Statistics: how often is damping active, what is the average COL sizing factor?
    m->hessDamped = 0;
    m->averageSizingFactor = 0.0;

    for (iBlock=0; iBlock<nBlocks; iBlock++) {
      nVarLocal = m->hess[iBlock].M();

      // smallGamma and smallDelta are subvectors of gamma and delta,
      // corresponding to partially separability
      smallGamma.Submatrix(m->gammaMat, nVarLocal, m->gammaMat.N(),
        blocks_[iBlock], 0);
      smallDelta.Submatrix(m->deltaMat, nVarLocal, m->deltaMat.N(),
        blocks_[iBlock], 0);

      // Is this the first iteration or the first after a Hessian reset?
      firstIter = (m->noUpdateCounter[iBlock] == -1);

      // Update sTs, sTs_ and sTy, sTy_
      m->deltaNormOld(iBlock) = m->deltaNorm(iBlock);
      m->deltaGammaOld(iBlock) = m->deltaGamma(iBlock);
      m->deltaNorm(iBlock) = adotb(smallDelta, smallDelta);
      m->deltaGamma(iBlock) = adotb(smallDelta, smallGamma);

      // Sizing before the update
      if (hessScaling < 4 && firstIter)
        sizeInitialHessian(m, smallGamma, smallDelta, iBlock, hessScaling);
      else if (hessScaling == 4)
        sizeHessianCOL(m, smallGamma, smallDelta, iBlock);

      // Compute the new update
      if (updateType == 1) {
        calcSR1(m, smallGamma, smallDelta, iBlock);

        // Prepare to compute fallback update as well
        m->hess = m->hess2;

        // Sizing the fallback update
        if (fallback_scaling_ < 4 && firstIter)
          sizeInitialHessian(m, smallGamma, smallDelta, iBlock, fallback_scaling_);
        else if (fallback_scaling_ == 4)
          sizeHessianCOL(m, smallGamma, smallDelta, iBlock);

        // Compute fallback update
        if (fallback_update_ == 2)
          calcBFGS(m, smallGamma, smallDelta, iBlock);

        // Reset pointer
        m->hess = m->hess1;
      } else if (updateType == 2) {
        calcBFGS(m, smallGamma, smallDelta, iBlock);
      }

      // If an update is skipped to often, reset Hessian block
      if (m->noUpdateCounter[iBlock] > max_consec_skipped_updates_) {
        resetHessian(m, iBlock);
      }
    }

    // statistics: average sizing factor
    m->averageSizingFactor /= nBlocks;
  }


  void Blocksqp::
  calcHessianUpdateLimitedMemory(BlocksqpMemory* m, int updateType, int hessScaling) const {
    int iBlock, nBlocks, nVarLocal;
    blocksqp::Matrix smallGamma, smallDelta;
    blocksqp::Matrix gammai, deltai;
    int i, m2, pos, posOldest, posNewest;
    int hessDamped, hessSkipped;
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

    for (iBlock=0; iBlock<nBlocks; iBlock++) {
      nVarLocal = m->hess[iBlock].M();

      // smallGamma and smallDelta are submatrices of gammaMat, deltaMat,
      // i.e. subvectors of gamma and delta from m prev. iterations
      smallGamma.Submatrix(m->gammaMat, nVarLocal, m->gammaMat.N(),
        blocks_[iBlock], 0);
      smallDelta.Submatrix(m->deltaMat, nVarLocal, m->deltaMat.N(),
        blocks_[iBlock], 0);

      // Memory structure
      if (m->itCount > smallGamma.N()) {
        m2 = smallGamma.N();
        posOldest = m->itCount % m2;
        posNewest = (m->itCount-1) % m2;
      } else {
        m2 = m->itCount;
        posOldest = 0;
        posNewest = m2-1;
      }

      // Set B_0 (pretend it's the first step)
      calcInitialHessian(m, iBlock);
      m->deltaNorm(iBlock) = 1.0;
      m->deltaNormOld(iBlock) = 1.0;
      m->deltaGamma(iBlock) = 0.0;
      m->deltaGammaOld(iBlock) = 0.0;
      m->noUpdateCounter[iBlock] = -1;

      // Size the initial update, but with the most recent delta/gamma-pair
      gammai.Submatrix(smallGamma, nVarLocal, 1, 0, posNewest);
      deltai.Submatrix(smallDelta, nVarLocal, 1, 0, posNewest);
      sizeInitialHessian(m, gammai, deltai, iBlock, hessScaling);

      for (i=0; i<m2; i++) {
        pos = (posOldest+i) % m2;

        // Get new vector from list
        gammai.Submatrix(smallGamma, nVarLocal, 1, 0, pos);
        deltai.Submatrix(smallDelta, nVarLocal, 1, 0, pos);

        // Update sTs, sTs_ and sTy, sTy_
        m->deltaNormOld(iBlock) = m->deltaNorm(iBlock);
        m->deltaGammaOld(iBlock) = m->deltaGamma(iBlock);
        m->deltaNorm(iBlock) = adotb(deltai, deltai);
        m->deltaGamma(iBlock) = adotb(gammai, deltai);

        // Save statistics, we want to record them only for the most recent update
        averageSizingFactor = m->averageSizingFactor;
        hessDamped = m->hessDamped;
        hessSkipped = m->hessSkipped;

        // Selective sizing before the update
        if (hessScaling == 4) sizeHessianCOL(m, gammai, deltai, iBlock);

        // Compute the new update
        if (updateType == 1) {
          calcSR1(m, gammai, deltai, iBlock);
        } else if (updateType == 2) {
          calcBFGS(m, gammai, deltai, iBlock);
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
      if (m->noUpdateCounter[iBlock] > max_consec_skipped_updates_) {
        resetHessian(m, iBlock);
      }
    }
    //blocks
    m->averageSizingFactor /= nBlocks;
  }


  void Blocksqp::
  calcBFGS(BlocksqpMemory* m, const blocksqp::Matrix &gamma,
    const blocksqp::Matrix &delta, int iBlock) const {
    int i, j, k, dim = gamma.M();
    blocksqp::Matrix Bdelta;
    blocksqp::SymMatrix *B;
    double h1 = 0.0;
    double h2 = 0.0;
    double thetaPowell = 0.0;
    int damped;

    /* Work with a local copy of gamma because damping may need to change gamma.
     * Note that m->gamma needs to remain unchanged!
     * This may be important in a limited memory context:
     * When information is "forgotten", B_i-1 is different and the
     *  original gamma might lead to an undamped update with the new B_i-1! */
    blocksqp::Matrix gamma2 = gamma;

    B = &m->hess[iBlock];

    // Bdelta = B*delta (if sizing is enabled, B is the sized B!)
    // h1 = delta^T * B * delta
    // h2 = delta^T * gamma
    Bdelta.Dimension(dim).Initialize(0.0);
    for (i=0; i<dim; i++) {
        for (k=0; k<dim; k++)
          Bdelta(i) += (*B)(i, k) * delta(k);

        h1 += delta(i) * Bdelta(i);
        //h2 += delta(i) * gamma(i);
      }
    h2 = m->deltaGamma(iBlock);

    /* Powell's damping strategy to maintain pos. def. (Nocedal/Wright p.537; SNOPT paper)
     * Interpolates between current approximation and unmodified BFGS */
    damped = 0;
    if (hess_damp_)
      if (h2 < hess_damp_fac_ * h1 / m->alpha && fabs(h1 - h2) > 1.0e-12) {
        // At the first iteration h1 and h2 are equal due to COL scaling

        thetaPowell = (1.0 - hess_damp_fac_)*h1 / (h1 - h2);

        // Redefine gamma and h2 = delta^T * gamma
        h2 = 0.0;
        for (i=0; i<dim; i++) {
          gamma2(i) = thetaPowell*gamma2(i) + (1.0 - thetaPowell)*Bdelta(i);
          h2 += delta(i) * gamma2(i);
        }

        // Also redefine deltaGamma for computation of sizing factor in the next iteration
        m->deltaGamma(iBlock) = h2;

        damped = 1;
      }

    // For statistics: count number of damped blocks
    m->hessDamped += damped;

    // B_k+1 = B_k - Bdelta * (Bdelta)^T / h1 + gamma * gamma^T / h2
    double myEps = 1.0e2 * eps_;
    if (fabs(h1) < myEps || fabs(h2) < myEps) {
      // don't perform update because of bad condition, might introduce negative eigenvalues
      m->noUpdateCounter[iBlock]++;
      m->hessDamped -= damped;
      m->hessSkipped++;
      m->nTotalSkippedUpdates++;
    } else {
      for (i=0; i<dim; i++)
        for (j=i; j<dim; j++)
          (*B)(i, j) = (*B)(i, j) - Bdelta(i) * Bdelta(j) / h1
            + gamma2(i) * gamma2(j) / h2;

      m->noUpdateCounter[iBlock] = 0;
    }
  }


  void Blocksqp::
  calcSR1(BlocksqpMemory* m, const blocksqp::Matrix &gamma, const blocksqp::Matrix &delta,
    int iBlock) const {
    int i, j, k, dim = gamma.M();
    blocksqp::Matrix gmBdelta;
    blocksqp::SymMatrix *B;
    double myEps = 1.0e2 * eps_;
    double r = 1.0e-8;
    double h = 0.0;

    B = &m->hess[iBlock];

    // gmBdelta = gamma - B*delta
    // h = (gamma - B*delta)^T * delta
    gmBdelta.Dimension(dim);
    for (i=0; i<dim; i++) {
      gmBdelta(i) = gamma(i);
      for (k=0; k<dim; k++)
        gmBdelta(i) -= ((*B)(i, k) * delta(k));

      h += (gmBdelta(i) * delta(i));
    }

    // B_k+1 = B_k + gmBdelta * gmBdelta^T / h
    if (fabs(h) < r * l2VectorNorm(delta) * l2VectorNorm(gmBdelta) || fabs(h) < myEps) {
      // Skip update if denominator is too small
      m->noUpdateCounter[iBlock]++;
      m->hessSkipped++;
      m->nTotalSkippedUpdates++;
    } else {
      for (i=0; i<dim; i++)
        for (j=i; j<dim; j++)
          (*B)(i, j) = (*B)(i, j) + gmBdelta(i) * gmBdelta(j) / h;
      m->noUpdateCounter[iBlock] = 0;
    }
  }


  /**
   * Set deltaXi and gamma as a column in the matrix containing
   * the m most recent delta and gamma
   */
  void Blocksqp::updateDeltaGamma(BlocksqpMemory* m) const {
    int nVar = m->gammaMat.M();
    int m2 = m->gammaMat.N();

    if (m2 == 1)
      return;

    m->deltaXi.Submatrix(m->deltaMat, nVar, 1, 0, m->itCount % m2);
    m->gamma.Submatrix(m->gammaMat, nVar, 1, 0, m->itCount % m2);
  }

  void Blocksqp::
  computeNextHessian(BlocksqpMemory* m, int idx, int maxQP) const {
    // Compute fallback update only once
    if (idx == 1) {
        // Switch storage
        m->hess = m->hess2;

        // If last block contains exact Hessian, we need to copy it
        if (which_second_derv_ == 1)
          for (int i=0; i<m->hess[nblocks_-1].M(); i++)
            for (int j=i; j<m->hess[nblocks_-1].N(); j++)
              m->hess2[nblocks_-1](i, j) = m->hess1[nblocks_-1](i, j);

        // Limited memory: compute fallback update only when needed
        if (hess_lim_mem_) {
            m->itCount--;
            int hessDampSave = hess_damp_;
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
        for (int iBlock=0; iBlock<nblocks_; iBlock++)
          for (int i=0; i<m->hess[iBlock].M(); i++)
            for (int j=i; j<m->hess[iBlock].N(); j++) {
                m->hess2[iBlock](i, j) *= mu;
                m->hess2[iBlock](i, j) += mu1 * m->hess1[iBlock](i, j);
              }
      }
  }


  /**
   * Inner loop of SQP algorithm:
   * Solve a sequence of QPs until pos. def. assumption (G3*) is satisfied.
   */
  int Blocksqp::
  solveQP(BlocksqpMemory* m, blocksqp::Matrix &deltaXi, blocksqp::Matrix &lambdaQP,
    bool matricesChanged) const {
    blocksqp::Matrix jacT;
    int maxQP, l;
    if (globalization_ == 1 &&
        hess_update_ == 1 &&
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
    qpOASES::Matrix *A;
    qpOASES::SymmetricMatrix *H;
    if (matricesChanged) {
        if (sparse_qp_) {
            A = new qpOASES::SparseMatrix(m->prob->nCon, m->prob->nVar,
                                           m->jacIndRow, m->jacIndCol, m->jacNz);
          } else {
            // transpose Jacobian (qpOASES needs row major arrays)
            Transpose(m->constrJac, jacT);
            A = new qpOASES::DenseMatrix(m->prob->nCon, m->prob->nVar, m->prob->nVar, jacT.ARRAY());
          }
      }
    double *g = m->gradObj.ARRAY();
    double *lb = m->deltaBl.ARRAY();
    double *lu = m->deltaBu.ARRAY();
    double *lbA = m->deltaBl.ARRAY() + m->prob->nVar;
    double *luA = m->deltaBu.ARRAY() + m->prob->nVar;

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

    if (maxQP > 1) {
        // Store last successful QP in temporary storage
        (*m->qpSave) = *m->qp;
        /** \todo Storing the active set would be enough but then the QP object
         *        must be properly reset after unsuccessful (SR1-)attempt.
         *        Moreover, passing a guessed active set doesn't yield
         *        exactly the same result as hotstarting a QP. This has
         *        something to do with how qpOASES handles user-given
         *        active sets (->check qpOASES source code). */
      }

    // Other variables for qpOASES
    double cpuTime = matricesChanged ? max_time_qp_ : 0.1*max_time_qp_;
    int maxIt = matricesChanged ? max_it_qp_ : 0.1*max_it_qp_;
    qpOASES::SolutionAnalysis solAna;
    qpOASES::returnValue ret;

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
          *m->qp = *m->qpSave;
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
            if (sparse_qp_) {
                // Convert block-Hessian to sparse format
                convertHessian(m, m->prob, eps_, m->hess, m->hessNz,
                                      m->hessIndRow, m->hessIndCol, m->hessIndLo);
                H = new qpOASES::SymSparseMat(m->prob->nVar, m->prob->nVar,
                                               m->hessIndRow, m->hessIndCol,
                                               m->hessNz);
                dynamic_cast<qpOASES::SymSparseMat*>(H)->createDiagInfo();
              } else {
                // Convert block-Hessian to double array
                convertHessian(m, m->prob, eps_, m->hess);
                H = new qpOASES::SymDenseMat(m->prob->nVar, m->prob->nVar,
                  m->prob->nVar, m->hessNz);
              }
          }

        /*
         * Call qpOASES
         */
        if (debug_level_ > 2) {
          dumpQPCpp(m, m->prob, m->qp, sparse_qp_);
        }
        if (matricesChanged) {
            maxIt = max_it_qp_;
            cpuTime = max_time_qp_;
            if (m->qp->getStatus() == qpOASES::QPS_HOMOTOPYQPSOLVED ||
                m->qp->getStatus() == qpOASES::QPS_SOLVED) {
                ret = m->qp->hotstart(H, g, A, lb, lu, lbA, luA, maxIt, &cpuTime);
              } else {
                ret = m->qp->init(H, g, A, lb, lu, lbA, luA, maxIt, &cpuTime);
              }
          } else if (!matricesChanged) {
            // Second order correction: H and A do not change
            maxIt = 0.1*max_it_qp_;
            cpuTime = 0.1*max_time_qp_;
            ret = m->qp->hotstart(g, lb, lu, lbA, luA, maxIt, &cpuTime);
          }

        /*
         * Check assumption (G3*) if nonconvex QP was solved
         */
        if (l < maxQP-1 && matricesChanged) {
            if (ret == qpOASES::SUCCESSFUL_RETURN) {
                if (sparse_qp_ == 2) {
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
    m->qp->getPrimalSolution(deltaXi.ARRAY());
    m->qp->getDualSolution(lambdaQP.ARRAY());
    m->qpObj = m->qp->getObjVal();

    // Compute constrJac*deltaXi, need this for second order correction step
    if (sparse_qp_) {
      Atimesb(m->jacNz, m->jacIndRow, m->jacIndCol, deltaXi, m->AdeltaXi);
    } else {
      Atimesb(m->constrJac, deltaXi, m->AdeltaXi);
    }

    // Print qpOASES error code, if any
    if (ret != qpOASES::SUCCESSFUL_RETURN && matricesChanged)
      printf("qpOASES error message: \"%s\"\n",
              qpOASES::getGlobalMessageHandler()->getErrorCodeMessage(ret));

    // Point Hessian again to the first Hessian
    m->hess = m->hess1;

    /* For full-memory Hessian: Restore fallback Hessian if convex combinations
     * were used during the loop */
    if (!hess_lim_mem_ && maxQP > 2 && matricesChanged) {
        double mu = 1.0 / l;
        double mu1 = 1.0 - mu;
        int nBlocks = (which_second_derv_ == 1) ? nblocks_-1 : nblocks_;
        for (int iBlock=0; iBlock<nBlocks; iBlock++)
          for (int i=0; i<m->hess[iBlock].M(); i++)
            for (int j=i; j<m->hess[iBlock].N(); j++) {
                m->hess2[iBlock](i, j) *= mu;
                m->hess2[iBlock](i, j) += mu1 * m->hess1[iBlock](i, j);
              }
      }

    /* Return code depending on qpOASES returnvalue
     * 0: Success
     * 1: Maximum number of iterations reached
     * 2: Unbounded
     * 3: Infeasible
     * 4: Other error */
    if (ret == qpOASES::SUCCESSFUL_RETURN)
      return 0;
    else if (ret == qpOASES::RET_MAX_NWSR_REACHED)
      return 1;
    else if (ret == qpOASES::RET_HESSIAN_NOT_SPD ||
             ret == qpOASES::RET_HESSIAN_INDEFINITE ||
             ret == qpOASES::RET_INIT_FAILED_UNBOUNDEDNESS ||
             ret == qpOASES::RET_QP_UNBOUNDED ||
             ret == qpOASES::RET_HOTSTART_STOPPED_UNBOUNDEDNESS)
      return 2;
    else if (ret == qpOASES::RET_INIT_FAILED_INFEASIBILITY ||
             ret == qpOASES::RET_QP_INFEASIBLE ||
             ret == qpOASES::RET_HOTSTART_STOPPED_INFEASIBILITY)
      return 3;
    else
      return 4;
  }


  /**
   * Set bounds on the step (in the QP), either according
   * to variable bounds in the NLP or according to
   * trust region box radius
   */
  void Blocksqp::updateStepBounds(BlocksqpMemory* m, bool soc) const {
    int i;
    int nVar = m->prob->nVar;
    int nCon = m->prob->nCon;

    // Bounds on step
    for (i=0; i<nVar; i++) {
      if (m->prob->bl(i) != inf)
        m->deltaBl(i) = m->prob->bl(i) - m->xi(i);
      else
        m->deltaBl(i) = inf;

      if (m->prob->bu(i) != inf)
        m->deltaBu(i) = m->prob->bu(i) - m->xi(i);
      else
        m->deltaBu(i) = inf;
    }

    // Bounds on linearized constraints
    for (i=0; i<nCon; i++) {
      if (m->prob->bl(nVar+i) != inf) {
        m->deltaBl(nVar+i) = m->prob->bl(nVar+i) - m->constr(i);
        if (soc) m->deltaBl(nVar+i) += m->AdeltaXi(i);
      } else {
        m->deltaBl(nVar+i) = inf;
      }

      if (m->prob->bu(nVar+i) != inf) {
        m->deltaBu(nVar+i) = m->prob->bu(nVar+i) - m->constr(i);
        if (soc) m->deltaBu(nVar+i) += m->AdeltaXi(i);
      } else {
        m->deltaBu(nVar+i) = inf;
      }
    }
  }

  void Blocksqp::
  printProgress(BlocksqpMemory* m, BlocksqpProblem *prob,
    bool hasConverged) const {
    /*
     * m->steptype:
     *-1: full step was accepted because it reduces the KKT error although line search failed
     * 0: standard line search step
     * 1: Hessian has been reset to identity
     * 2: feasibility restoration heuristic has been called
     * 3: feasibility restoration phase has been called
     */

    if (m->itCount == 0) {
      if (print_level_ > 0) {
        prob->printInfo();

        // Headline
        printf("%-8s", "   it");
        printf("%-21s", " qpIt");
        printf("%-9s", "obj");
        printf("%-11s", "feas");
        printf("%-7s", "opt");
        if (print_level_ > 1) {
          printf("%-11s", "|lgrd|");
          printf("%-9s", "|stp|");
          printf("%-10s", "|lstp|");
        }
        printf("%-8s", "alpha");
        if (print_level_ > 1) {
          printf("%-6s", "nSOCS");
          printf("%-18s", "sk, da, sca");
          printf("%-6s", "QPr,mu");
        }
        printf("\n");

        // Values for first iteration
        printf("%5i  ", m->itCount);
        printf("%11i ", 0);
        printf("% 10e  ", m->obj);
        printf("%-10.2e", m->cNormS);
        printf("%-10.2e", m->tol);
        printf("\n");
      }

      if (debug_level_ > 0) {
        // Print everything in a CSV file as well
        fprintf(m->progressFile, "%23.16e, %23.16e, %23.16e, %23.16e, %23.16e, %23.16e, "
                "%23.16e, %23.16e, %i, %i, %23.16e, %i, %23.16e\n",
                 m->obj, m->cNormS, m->tol, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0.0, 0, 0.0);
      }
    } else {
      // Every twenty iterations print headline
      if (m->itCount % 20 == 0 && print_level_ > 0) {
        printf("%-8s", "   it");
        printf("%-21s", " qpIt");
        printf("%-9s", "obj");
        printf("%-11s", "feas");
        printf("%-7s", "opt");
        if (print_level_ > 1) {
            printf("%-11s", "|lgrd|");
            printf("%-9s", "|stp|");
            printf("%-10s", "|lstp|");
          }
        printf("%-8s", "alpha");
        if (print_level_ > 1) {
            printf("%-6s", "nSOCS");
            printf("%-18s", "sk, da, sca");
            printf("%-6s", "QPr,mu");
          }
        printf("\n");
      }

      // All values
      if (print_level_ > 0) {
        printf("%5i  ", m->itCount);
        printf("%5i+%5i ", m->qpIterations, m->qpIterations2);
        printf("% 10e  ", m->obj);
        printf("%-10.2e", m->cNormS);
        printf("%-10.2e", m->tol);
        if (print_level_ > 1) {
            printf("%-10.2e", m->gradNorm);
            printf("%-10.2e", lInfVectorNorm(m->deltaXi));
            printf("%-10.2e", m->lambdaStepNorm);
        }
        printf("%-9.1e", m->alpha);

        if (print_level_ > 1) {
          printf("%5i", m->nSOCS);
          printf("%3i, %3i, %-9.1e", m->hessSkipped, m->hessDamped, m->averageSizingFactor);
          printf("%i, %-9.1e", m->qpResolve, l1VectorNorm(m->deltaH)/nblocks_);
        }
        printf("\n");
      }

      if (debug_level_ > 0) {
        // Print everything in a CSV file as well
        fprintf(m->progressFile, "%23.16e, %23.16e, %23.16e, %23.16e, %23.16e, "
                "%23.16e, %23.16e, %i, %i, %i, %23.16e, %i, %23.16e\n",
                 m->obj, m->cNormS, m->tol, m->gradNorm,
                 lInfVectorNorm(m->deltaXi),
                 m->lambdaStepNorm, m->alpha, m->nSOCS, m->hessSkipped,
                 m->hessDamped, m->averageSizingFactor,
                 m->qpResolve, l1VectorNorm(m->deltaH)/nblocks_);

        // Print update sequence
        fprintf(m->updateFile, "%i\t", m->qpResolve);
      }
    }

    // Print Debug information
    printDebug(m);

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

    if (print_level_ > 0) {
      if (hasConverged && m->steptype < 2) {
        printf("\n***CONVERGENCE ACHIEVED!***\n");
      }
    }
  }


  void Blocksqp::initStats(BlocksqpMemory* m) const {
    blocksqp::PATHSTR filename;

    // Open files

    if (debug_level_ > 0) {
      // SQP progress
      strcpy(filename, m->outpath);
      strcat(filename, "sqpits.csv");
      m->progressFile = fopen(filename, "w");

      // Update sequence
      strcpy(filename, m->outpath);
      strcat(filename, "updatesequence.txt");
      m->updateFile = fopen(filename, "w");
    }

    if (debug_level_ > 1) {
      // Primal variables
      strcpy(filename, m->outpath);
      strcat(filename, "pv.csv");
      m->primalVarsFile = fopen(filename, "w");

      // Dual variables
      strcpy(filename, m->outpath);
      strcat(filename, "dv.csv");
      m->dualVarsFile = fopen(filename, "w");
    }

    m->itCount = 0;
    m->qpItTotal = 0;
    m->qpIterations = 0;
    m->hessSkipped = 0;
    m->hessDamped = 0;
    m->averageSizingFactor = 0.0;
  }


  void Blocksqp::printPrimalVars(BlocksqpMemory* m, const blocksqp::Matrix &xi) const {
    for (int i=0; i<xi.M()-1; i++)
      fprintf(m->primalVarsFile, "%23.16e ", xi(i));
    fprintf(m->primalVarsFile, "%23.16e\n", xi(xi.M()-1));
  }


  void Blocksqp::printDualVars(BlocksqpMemory* m, const blocksqp::Matrix &lambda) const {
    for (int i=0; i<lambda.M()-1; i++)
      fprintf(m->dualVarsFile, "%23.16e ", lambda(i));
    fprintf(m->dualVarsFile, "%23.16e\n", lambda(lambda.M()-1));
  }


  void Blocksqp::
  printHessian(BlocksqpMemory* m, int nBlocks, blocksqp::SymMatrix *&hess) const {
    blocksqp::PATHSTR filename;
    int offset, i, j, iBlock, nVar;

    nVar = 0;
    for (iBlock=0; iBlock<nBlocks; iBlock++)
      nVar += hess[iBlock].M();

    blocksqp::SymMatrix fullHessian;
    fullHessian.Dimension(nVar).Initialize(0.0);

    strcpy(filename, m->outpath);
    strcat(filename, "hes.m");
    m->hessFile = fopen(filename, "w");

    offset = 0;
    for (iBlock=0; iBlock<nBlocks; iBlock++) {
        for (i=0; i<hess[iBlock].N(); i++)
          for (j=i; j<hess[iBlock].N(); j++)
            fullHessian(offset + i, offset + j) = hess[iBlock](i, j);

        offset += hess[iBlock].N();
      }

    fprintf(m->hessFile, "H=");
    fullHessian.Print(m->hessFile, 23, 1);
    fprintf(m->hessFile, "\n");
    fclose(m->hessFile);
  }


  void Blocksqp::
  printHessian(BlocksqpMemory* m, int nVar, double *hesNz, int *hesIndRow, int *hesIndCol) const {
    blocksqp::PATHSTR filename;

    strcpy(filename, m->outpath);
    strcat(filename, "hes.dat");
    m->hessFile = fopen(filename, "w");

    printSparseMatlab(m, m->hessFile, nVar, nVar, hesNz, hesIndRow, hesIndCol);

    fprintf(m->hessFile, "\n");
    fclose(m->hessFile);
  }


  void Blocksqp::printJacobian(BlocksqpMemory* m, const blocksqp::Matrix &constrJac) const {
    blocksqp::PATHSTR filename;

    strcpy(filename, m->outpath);
    strcat(filename, "jac.m");
    m->jacFile = fopen(filename, "w");

    fprintf(m->jacFile, "A=");
    constrJac.Print(m->jacFile, 23, 1);
    fprintf(m->jacFile, "\n");

    fclose(m->jacFile);
  }


  void Blocksqp::
  printJacobian(BlocksqpMemory* m, int nCon, int nVar, double *jacNz,
                int *jacIndRow, int *jacIndCol) const {
    blocksqp::PATHSTR filename;

    strcpy(filename, m->outpath);
    strcat(filename, "jac.dat");
    m->jacFile = fopen(filename, "w");

    printSparseMatlab(m, m->jacFile, nCon, nVar, jacNz, jacIndRow, jacIndCol);

    fprintf(m->jacFile, "\n");
    fclose(m->jacFile);
  }

  void Blocksqp::
  printSparseMatlab(BlocksqpMemory* m, FILE *file, int nRow, int nCol,
                    double *nz, int *indRow, int *indCol) const {
    int i, j, count;

    count = 0;
    fprintf(file, "%i %i 0\n", nRow, nCol);
    for (i=0; i<nCol; i++)
      for (j=indCol[i]; j<indCol[i+1]; j++) {
          // +1 for MATLAB indices!
          fprintf(file, "%i %i %23.16e\n", indRow[count]+1, i+1, nz[count]);
          count++;
        }
  }


  void Blocksqp::
  printDebug(BlocksqpMemory* m) const {
    if (debug_level_ > 1) {
        printPrimalVars(m, m->xi);
        printDualVars(m, m->lambda);
      }
  }

  void Blocksqp::printCppNull(BlocksqpMemory* m, FILE *outfile, char* varname) const {
    fprintf(outfile, "    double *%s = 0;\n", varname);
  }


  void Blocksqp::
  printVectorCpp(BlocksqpMemory* m, FILE *outfile, double *vec, int len, char* varname) const {
    int i;

    fprintf(outfile, "    double %s[%i] = { ", varname, len);
    for (i=0; i<len; i++) {
      fprintf(outfile, "%23.16e", vec[i]);
      if (i != len-1)
        fprintf(outfile, ", ");
      if ((i+1) % 10 == 0)
        fprintf(outfile, "\n          ");
    }
    fprintf(outfile, " };\n\n");
  }


  void Blocksqp::printVectorCpp(BlocksqpMemory* m, FILE *outfile, int *vec,
    int len, char* varname) const {
    int i;

    fprintf(outfile, "    int %s[%i] = { ", varname, len);
    for (i=0; i<len; i++) {
      fprintf(outfile, "%i", vec[i]);
      if (i != len-1)
        fprintf(outfile, ", ");
      if ((i+1) % 15 == 0)
        fprintf(outfile, "\n          ");
    }
    fprintf(outfile, " };\n\n");
  }


  void Blocksqp::
  dumpQPCpp(BlocksqpMemory* m, BlocksqpProblem *prob,
    qpOASES::SQProblem *qp, int sparseQP) const {
    int i, j;
    blocksqp::PATHSTR filename;
    FILE *outfile;
    int n = prob->nVar;

    // Print dimensions
    strcpy(filename, m->outpath);
    strcat(filename, "qpoases_dim.dat");
    outfile = fopen(filename, "w");
    fprintf(outfile, "%i %i\n", n, prob->nCon);
    fclose(outfile);

    // Print Hessian
    if (sparseQP) {
      strcpy(filename, m->outpath);
      strcat(filename, "qpoases_H_sparse.dat");
      outfile = fopen(filename, "w");
      for (i=0; i<prob->nVar+1; i++)
        fprintf(outfile, "%i ", m->hessIndCol[i]);
      fprintf(outfile, "\n");

      for (i=0; i<m->hessIndCol[prob->nVar]; i++)
        fprintf(outfile, "%i ", m->hessIndRow[i]);
      fprintf(outfile, "\n");

      for (i=0; i<m->hessIndCol[prob->nVar]; i++)
        fprintf(outfile, "%23.16e ", m->hessNz[i]);
      fprintf(outfile, "\n");
      fclose(outfile);
    }
    strcpy(filename, m->outpath);
    strcat(filename, "qpoases_H.dat");
    outfile = fopen(filename, "w");
    int blockCnt = 0;
    for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
        if (i == blocks_[blockCnt+1]) blockCnt++;
        if (j >= blocks_[blockCnt] && j < blocks_[blockCnt+1]) {
          fprintf(outfile, "%23.16e ", m->hess[blockCnt](i - blocks_[blockCnt],
            j - blocks_[blockCnt]));
        } else {
          fprintf(outfile, "0.0 ");
        }
      }
      fprintf(outfile, "\n");
    }
    fclose(outfile);

    // Print gradient
    strcpy(filename, m->outpath);
    strcat(filename, "qpoases_g.dat");
    outfile = fopen(filename, "w");
    for (i=0; i<n; i++)
      fprintf(outfile, "%23.16e ", m->gradObj(i));
    fprintf(outfile, "\n");
    fclose(outfile);

    // Print Jacobian
    strcpy(filename, m->outpath);
    strcat(filename, "qpoases_A.dat");
    outfile = fopen(filename, "w");
    if (sparseQP) {
      // Always print dense Jacobian
      blocksqp::Matrix constrJacTemp;
      constrJacTemp.Dimension(prob->nCon, prob->nVar).Initialize(0.0);
      for (i=0; i<prob->nVar; i++)
        for (j=m->jacIndCol[i]; j<m->jacIndCol[i+1]; j++)
          constrJacTemp(m->jacIndRow[j], i) = m->jacNz[j];
      for (i=0; i<prob->nCon; i++) {
        for (j=0; j<n; j++)
          fprintf(outfile, "%23.16e ", constrJacTemp(i, j));
        fprintf(outfile, "\n");
      }
      fclose(outfile);
    } else {
      for (i=0; i<prob->nCon; i++) {
        for (j=0; j<n; j++)
          fprintf(outfile, "%23.16e ", m->constrJac(i, j));
        fprintf(outfile, "\n");
      }
      fclose(outfile);
    }

    if (sparseQP) {
      strcpy(filename, m->outpath);
      strcat(filename, "qpoases_A_sparse.dat");
      outfile = fopen(filename, "w");
      for (i=0; i<prob->nVar+1; i++)
        fprintf(outfile, "%i ", m->jacIndCol[i]);
      fprintf(outfile, "\n");

      for (i=0; i<m->jacIndCol[prob->nVar]; i++)
        fprintf(outfile, "%i ", m->jacIndRow[i]);
      fprintf(outfile, "\n");

      for (i=0; i<m->jacIndCol[prob->nVar]; i++)
        fprintf(outfile, "%23.16e ", m->jacNz[i]);
      fprintf(outfile, "\n");
      fclose(outfile);
    }

    // Print variable lower bounds
    strcpy(filename, m->outpath);
    strcat(filename, "qpoases_lb.dat");
    outfile = fopen(filename, "w");
    for (i=0; i<n; i++)
      fprintf(outfile, "%23.16e ", m->deltaBl(i));
    fprintf(outfile, "\n");
    fclose(outfile);

    // Print variable upper bounds
    strcpy(filename, m->outpath);
    strcat(filename, "qpoases_ub.dat");
    outfile = fopen(filename, "w");
    for (i=0; i<n; i++)
      fprintf(outfile, "%23.16e ", m->deltaBu(i));
    fprintf(outfile, "\n");
    fclose(outfile);

    // Print constraint lower bounds
    strcpy(filename, m->outpath);
    strcat(filename, "qpoases_lbA.dat");
    outfile = fopen(filename, "w");
    for (i=0; i<prob->nCon; i++)
      fprintf(outfile, "%23.16e ", m->deltaBl(i+n));
    fprintf(outfile, "\n");
    fclose(outfile);

    // Print constraint upper bounds
    strcpy(filename, m->outpath);
    strcat(filename, "qpoases_ubA.dat");
    outfile = fopen(filename, "w");
    for (i=0; i<prob->nCon; i++)
      fprintf(outfile, "%23.16e ", m->deltaBu(i+n));
    fprintf(outfile, "\n");
    fclose(outfile);

    // Print active set
    qpOASES::Bounds b;
    qpOASES::Constraints c;
    qp->getBounds(b);
    qp->getConstraints(c);

    strcpy(filename, m->outpath);
    strcat(filename, "qpoases_as.dat");
    outfile = fopen(filename, "w");
    for (i=0; i<n; i++)
      fprintf(outfile, "%i ", b.getStatus(i));
    fprintf(outfile, "\n");
    for (i=0; i<prob->nCon; i++)
      fprintf(outfile, "%i ", c.getStatus(i));
    fprintf(outfile, "\n");
    fclose(outfile);
  }

  void Blocksqp::dumpQPMatlab(BlocksqpMemory* m, BlocksqpProblem *prob,
    int sparseQP) const {
    blocksqp::Matrix temp;
    blocksqp::PATHSTR filename;
    FILE *qpFile;
    FILE *vecFile;

    // Print vectors g, lb, lu, lbA, luA
    strcpy(filename, m->outpath);
    strcat(filename, "vec.m");
    vecFile = fopen(filename, "w");

    fprintf(vecFile, "g=");
    m->gradObj.Print(vecFile, 23, 1);
    fprintf(vecFile, "\n\n");

    temp.Submatrix(m->deltaBl, prob->nVar, 1, 0, 0);
    fprintf(vecFile, "lb=");
    temp.Print(vecFile, 23, 1);
    fprintf(vecFile, "\n\n");

    temp.Submatrix(m->deltaBu, prob->nVar, 1, 0, 0);
    fprintf(vecFile, "lu=");
    temp.Print(vecFile, 23, 1);
    fprintf(vecFile, "\n\n");

    temp.Submatrix(m->deltaBl, prob->nCon, 1, prob->nVar, 0);
    fprintf(vecFile, "lbA=");
    temp.Print(vecFile, 23, 1);
    fprintf(vecFile, "\n\n");

    temp.Submatrix(m->deltaBu, prob->nCon, 1, prob->nVar, 0);
    fprintf(vecFile, "luA=");
    temp.Print(vecFile, 23, 1);
    fprintf(vecFile, "\n");

    fclose(vecFile);

    // Print sparse Jacobian and Hessian
    if (sparseQP) {
        printJacobian(m, prob->nCon, prob->nVar, m->jacNz, m->jacIndRow, m->jacIndCol);
        printHessian(m, prob->nVar, m->hessNz, m->hessIndRow, m->hessIndCol);
      }

    // Print a script that correctly reads everything
    strcpy(filename, m->outpath);
    strcat(filename, "getqp.m");
    qpFile = fopen(filename, "w");

    fprintf(qpFile, "%% Read vectors g, lb, lu, lbA, luA\n");
    fprintf(qpFile, "vec;\n");
    fprintf(qpFile, "%% Read sparse Jacobian\n");
    fprintf(qpFile, "load jac.dat\n");
    fprintf(qpFile, "if jac(1) == 0\n");
    fprintf(qpFile, "    A = [];\n");
    fprintf(qpFile, "else\n");
    fprintf(qpFile, "    A = spconvert(jac);\n");
    fprintf(qpFile, "end\n");
    fprintf(qpFile, "%% Read sparse Hessian\n");
    fprintf(qpFile, "load hes.dat\n");
    fprintf(qpFile, "H = spconvert(hes);\n");

    fclose(qpFile);
  }

  /**
   * Allocate memory for variables
   * required by all optimization
   * algorithms except for the Jacobian
   */
  void Blocksqp::allocMin(BlocksqpMemory* m, BlocksqpProblem *prob) const {
    // current iterate
    m->xi.Dimension(prob->nVar).Initialize(0.0);

    // dual variables (for general constraints and variable bounds)
    m->lambda.Dimension(prob->nVar + prob->nCon).Initialize(0.0);

    // constraint vector with lower and upper bounds
    // (Box constraints are not included in the constraint list)
    m->constr.Dimension(prob->nCon).Initialize(0.0);

    // gradient of objective
    m->gradObj.Dimension(prob->nVar).Initialize(0.0);

    // gradient of Lagrangian
    m->gradLagrange.Dimension(prob->nVar).Initialize(0.0);
  }


  void Blocksqp::allocHess(BlocksqpMemory* m) const {
    int iBlock, varDim;

    // Create one Matrix for one diagonal block in the Hessian
    m->hess1 = new blocksqp::SymMatrix[nblocks_];
    for (iBlock=0; iBlock<nblocks_; iBlock++) {
        varDim = blocks_[iBlock+1] - blocks_[iBlock];
        m->hess1[iBlock].Dimension(varDim).Initialize(0.0);
      }

    // For SR1 or finite differences, maintain two Hessians
    if (hess_update_ == 1 || hess_update_ == 4) {
      m->hess2 = new blocksqp::SymMatrix[nblocks_];
      for (iBlock=0; iBlock<nblocks_; iBlock++) {
        varDim = blocks_[iBlock+1] - blocks_[iBlock];
        m->hess2[iBlock].Dimension(varDim).Initialize(0.0);
      }
    }

    // Set Hessian pointer
    m->hess = m->hess1;
  }

  /**
   * Convert diagonal block Hessian to double array.
   * Assumes that hessNz is already allocated.
   */
  void Blocksqp::convertHessian(BlocksqpMemory* m, BlocksqpProblem *prob,
    double eps, blocksqp::SymMatrix *&hess_) const {
    if (m->hessNz == 0) return;
    int count = 0;
    int blockCnt = 0;
    for (int i=0; i<prob->nVar; i++)
      for (int j=0; j<prob->nVar; j++) {
          if (i == blocks_[blockCnt+1])
            blockCnt++;
          if (j >= blocks_[blockCnt] && j < blocks_[blockCnt+1])
            m->hessNz[count++] = m->hess[blockCnt](i - blocks_[blockCnt],
              j - blocks_[blockCnt]);
          else
            m->hessNz[count++] = 0.0;
        }
  }

  /**
   * Convert array *hess to a single symmetric sparse matrix in
   * Harwell-Boeing format (as used by qpOASES)
   */
  void Blocksqp::
  convertHessian(BlocksqpMemory* m, BlocksqpProblem *prob, double eps,
                 blocksqp::SymMatrix *&hess_, double *&hessNz_,
                 int *&hessIndRow_, int *&hessIndCol_, int *&hessIndLo_) const {
    int iBlock, count, colCountTotal, rowOffset, i, j;
    int nnz, nCols, nRows;

    // 1) count nonzero elements
    nnz = 0;
    for (iBlock=0; iBlock<nblocks_; iBlock++)
      for (i=0; i<hess_[iBlock].N(); i++)
        for (j=i; j<hess_[iBlock].N(); j++)
          if (fabs(hess_[iBlock](i, j)) > eps) {
            nnz++;
            if (i != j) {
              // off-diagonal elements count twice
              nnz++;
            }
          }

    if (hessNz_ != 0) delete[] hessNz_;
    if (hessIndRow_ != 0) delete[] hessIndRow_;

    hessNz_ = new double[nnz];
    hessIndRow_ = new int[nnz + (prob->nVar+1) + prob->nVar];
    hessIndCol_ = hessIndRow_ + nnz;
    hessIndLo_ = hessIndCol_ + (prob->nVar+1);

    // 2) store matrix entries columnwise in hessNz
    count = 0; // runs over all nonzero elements
    colCountTotal = 0; // keep track of position in large matrix
    rowOffset = 0;
    for (iBlock=0; iBlock<nblocks_; iBlock++) {
      nCols = hess_[iBlock].N();
      nRows = hess_[iBlock].M();

      for (i=0; i<nCols; i++) {
        // column 'colCountTotal' starts at element 'count'
        hessIndCol_[colCountTotal] = count;

        for (j=0; j<nRows; j++)
          if (fabs(hess_[iBlock](i, j)) > eps) {
              hessNz_[count] = hess_[iBlock](i, j);
              hessIndRow_[count] = j + rowOffset;
              count++;
            }
        colCountTotal++;
      }

      rowOffset += nRows;
    }
    hessIndCol_[colCountTotal] = count;

    // 3) Set reference to lower triangular matrix
    for (j=0; j<prob->nVar; j++) {
      for (i=hessIndCol_[j]; i<hessIndCol_[j+1] && hessIndRow_[i]<j; i++) {}
      hessIndLo_[j] = i;
    }

    if (count != nnz)
      printf("Error in convertHessian: %i elements processed, "
            "should be %i elements!\n", count, nnz);
  }


  /**
   * Allocate memory for additional variables
   * needed by the algorithm
   */
  void Blocksqp::
  allocAlg(BlocksqpMemory* m, BlocksqpProblem *prob) const {
    int iBlock;
    int nVar = prob->nVar;
    int nCon = prob->nCon;

    // current step
    m->deltaMat.Dimension(nVar, hess_memsize_, nVar).Initialize(0.0);
    m->deltaXi.Submatrix(m->deltaMat, nVar, 1, 0, 0);
    // trial step (temporary variable, for line search)
    m->trialXi.Dimension(nVar, 1, nVar).Initialize(0.0);

    // bounds for step (QP subproblem)
    m->deltaBl.Dimension(nVar+nCon).Initialize(0.0);
    m->deltaBu.Dimension(nVar+nCon).Initialize(0.0);

    // product of constraint Jacobian with step (deltaXi)
    m->AdeltaXi.Dimension(nCon).Initialize(0.0);

    // dual variables of QP (simple bounds and general constraints)
    m->lambdaQP.Dimension(nVar+nCon).Initialize(0.0);

    // line search parameters
    m->deltaH.Dimension(nblocks_).Initialize(0.0);

    // filter as a set of pairs
    m->filter = new std::set< std::pair<double, double> >;

    // difference of Lagrangian gradients
    m->gammaMat.Dimension(nVar, hess_memsize_, nVar).Initialize(0.0);
    m->gamma.Submatrix(m->gammaMat, nVar, 1, 0, 0);

    // Scalars that are used in various Hessian update procedures
    m->noUpdateCounter = new int[nblocks_];
    for (iBlock=0; iBlock<nblocks_; iBlock++)
      m->noUpdateCounter[iBlock] = -1;

    // For selective sizing: for each block save sTs, sTs_, sTy, sTy_
    m->deltaNorm.Dimension(nblocks_).Initialize(1.0);
    m->deltaNormOld.Dimension(nblocks_).Initialize(1.0);
    m->deltaGamma.Dimension(nblocks_).Initialize(0.0);
    m->deltaGammaOld.Dimension(nblocks_).Initialize(0.0);
  }


  void Blocksqp::
  initIterate(BlocksqpMemory* m) const {
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

} // namespace casadi
