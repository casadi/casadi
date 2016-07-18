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

  BlocksqpProblem::BlocksqpProblem(const BlocksqpInterface& self, BlocksqpMemory* m)
    : self(self), m(m) {
    nVar = self.nx_;
    nCon = self.ng_;
    nBlocks = self.blocks_.size()-1;
    blockIdx = new int[nBlocks+1];
    for (int i=0; i<nBlocks+1; i++) blockIdx[i] = self.blocks_[i];

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

  extern "C"
  int CASADI_NLPSOL_BLOCKSQP_EXPORT
  casadi_register_nlpsol_blocksqp(Nlpsol::Plugin* plugin) {
    plugin->creator = BlocksqpInterface::creator;
    plugin->name = "blocksqp";
    plugin->doc = BlocksqpInterface::meta_doc.c_str();
    plugin->version = 30;
    return 0;
  }

  extern "C"
  void CASADI_NLPSOL_BLOCKSQP_EXPORT casadi_load_nlpsol_blocksqp() {
    Nlpsol::registerPlugin(casadi_register_nlpsol_blocksqp);
  }

  BlocksqpInterface::BlocksqpInterface(const std::string& name, const Function& nlp)
    : Nlpsol(name, nlp) {
  }


  BlocksqpInterface::~BlocksqpInterface() {
    clear_memory();
  }

  Options BlocksqpInterface::options_
  = {{&Nlpsol::options_},
     {{"blocksqp",
       {OT_DICT,
        "Options to be passed to BLOCKSQP"}},
      {"contype",
       {OT_INTVECTOR,
        "Type of constraint"}}
     }
  };

  void BlocksqpInterface::init(const Dict& opts) {
    // Call the init method of the base class
    Nlpsol::init(opts);

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="blocksqp") {
        //opts_ = op.second;
      } else if (op.first=="contype") {
        //contype_ = op.second;
      }
    }

    // Setup NLP functions
    create_function("nlp_fg", {"x", "p"}, {"f", "g"});
    Function gf_jg = create_function("nlp_gf_jg", {"x", "p"},
                                    {"f", "g", "grad:f:x", "jac:g:x"});
    sp_jac_ = gf_jg.sparsity_out("jac_g_x");

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

    // Allocate memory
    alloc_w(sp_jac_.nnz(), true); // jac
  }

  void BlocksqpInterface::init_memory(void* mem) const {
    Nlpsol::init_memory(mem);
    auto m = static_cast<BlocksqpMemory*>(mem);
  }

  void BlocksqpInterface::set_work(void* mem, const double**& arg, double**& res,
                                   int*& iw, double*& w) const {
    auto m = static_cast<BlocksqpMemory*>(mem);

    // Set work in base classes
    Nlpsol::set_work(mem, arg, res, iw, w);

    // Temporary memory
    m->jac = w; w += sp_jac_.nnz();
  }

  void BlocksqpInterface::solve(void* mem) const {
    auto m = static_cast<BlocksqpMemory*>(mem);

    int ret = 0;
    BlocksqpProblem *prob;
    blocksqp::SQPMethod *meth;
    blocksqp::SQPoptions *opts;
    blocksqp::SQPstats *stats;
    char outpath[255];
    strcpy(outpath, "./");

    // Create problem evaluation object
    vector<int> blocks = blocks_;
    prob = new BlocksqpProblem(*this, m);

    /*------------------------*/
    /* Options for SQP solver */
    /*------------------------*/
    opts = new blocksqp::SQPoptions();
    opts->opttol = 1.0e-12;
    opts->nlinfeastol = 1.0e-12;

    // 0: no globalization, 1: filter line search
    opts->globalization = 0;
    // 0: (scaled) identity, 1: SR1, 2: BFGS
    opts->hessUpdate = 0;
    // 0: initial Hessian is diagonal matrix, 1: scale initial Hessian according to Nocedal p.143,
    // 2: scale initial Hessian with Oren-Luenberger factor 3: scale initial Hessian
    //    with geometric mean of 1 and 2
    // 4: scale Hessian in every step with centered Oren-Luenberger sizing according to Tapia paper
    opts->hessScaling = 0;
    // scaling strategy for fallback BFGS update if SR1 and globalization is used
    opts->fallbackScaling = 0;
    // Size of limited memory
    opts->hessLimMem = 0;
    // If too many updates are skipped, reset Hessian
    opts->maxConsecSkippedUpdates = 200;
    // 0: full space Hessian approximation (ignore block structure), 1: blockwise updates
    opts->blockHess = 0;
    opts->whichSecondDerv = 0;
    opts->sparseQP = 1;
    opts->printLevel = 2;


    /*-------------------------------------------------*/
    /* Create blockSQP method object and run algorithm */
    /*-------------------------------------------------*/
    stats = new blocksqp::SQPstats( outpath );
    meth = new blocksqp::SQPMethod( prob, opts, stats );

    meth->init();
    ret = meth->run(100);
    meth->finish();
    if (ret==1) casadi_warning("Maximum number of iterations reached");

    // Get ptimal cost
    if (m->f) *m->f = meth->vars->obj;
    // Get primal solution
    casadi_copy(meth->vars->xi.array, nx_, m->x);
    // Get dual solution (simple bounds)
    if (m->lam_x) {
      casadi_copy(meth->vars->lambda.array, nx_, m->lam_x);
      casadi_scal(nx_, -1., m->lam_x);
    }
    // Get dual solution (nonlinear bounds)
    if (m->lam_g) {
      casadi_copy(meth->vars->lambda.array + nx_, ng_, m->lam_g);
      casadi_scal(ng_, -1., m->lam_g);
    }

    // Clean up
    delete prob;
    delete stats;
    delete opts;
    delete meth;
  }

} // namespace casadi
