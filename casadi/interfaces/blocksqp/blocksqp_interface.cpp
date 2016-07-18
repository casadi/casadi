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


  void BlocksqpProblem::convertJacobian(const blocksqp::Matrix &constrJac, double *&jacNz,
                                  int *&jacIndRow, int *&jacIndCol, bool firstCall) {
    int nnz, count, i, j;

    if (firstCall) {
      // 1st run: count nonzeros
      nnz = 0;
      for (j=0; j<nVar; j++)
        for (i=0; i<nCon; i++)
          if (fabs(constrJac(i, j) < inf)) nnz++;

      if (jacNz != 0) delete[] jacNz;
      if (jacIndRow != 0) delete[] jacIndRow;

      jacNz = new double[nnz];
      jacIndRow = new int[nnz + (nVar+1)];
      jacIndCol = jacIndRow + nnz;
    } else {
      nnz = jacIndCol[nVar];
      /* arrays jacInd* are already allocated! */
    }

    // 2nd run: store matrix entries columnwise in jacNz
    count = 0;
    for (j=0; j<nVar; j++) {
      jacIndCol[j] = count;
      for (i=0; i<nCon; i++)
        if (fabs(constrJac(i, j) < inf)) {
          jacNz[count] = constrJac(i, j);
          jacIndRow[count] = i;
          count++;
        }
    }
    jacIndCol[nVar] = count;
    if (count != nnz)
      printf("Error in convertJacobian: %i elements processed, "
             "should be %i elements!\n", count, nnz);
  }


  void BlocksqpProblem::initialize(blocksqp::Matrix &xi, blocksqp::Matrix &lambda,
                             blocksqp::Matrix &constrJac) {
    // set initial values for xi and lambda
    lambda.Initialize(0.0);
    for (int i=0; i<nVar; i++)
      xi(i) = m->x0 ? m->x0[i] : 0;
  }


  void BlocksqpProblem::initialize(blocksqp::Matrix &xi, blocksqp::Matrix &lambda,
                             double *&jacNz, int *&jacIndRow, int *&jacIndCol) {
    blocksqp::Matrix constrDummy, gradObjDummy, constrJac;
    blocksqp::SymMatrix *hessDummy;
    double objvalDummy;
    int info;

    // set initial values for xi and lambda
    lambda.Initialize(0.0);
    for (int i=0; i<nVar; i++)
      xi(i) = m->x0 ? m->x0[i] : 0;

    // find out Jacobian sparsity pattern by evaluating derivatives once
    constrDummy.Dimension(nCon).Initialize(0.0);
    gradObjDummy.Dimension(nVar).Initialize(0.0);
    constrJac.Dimension(nCon, nVar).Initialize(inf);
    evaluate(xi, lambda, &objvalDummy, constrDummy, gradObjDummy, constrJac,
             hessDummy, 1, &info);

    // allocate sparse Jacobian structures
    convertJacobian(constrJac, jacNz, jacIndRow, jacIndCol, 1);
  }

  /*
   * PROBLEM-SPECIFIC PART STARTS HERE
   */

  void BlocksqpProblem::evaluate(const blocksqp::Matrix &xi, const blocksqp::Matrix &lambda,
                           double *objval, blocksqp::Matrix &constr,
                           blocksqp::Matrix &gradObj, double *&jacNz, int *&jacIndRow,
                           int *&jacIndCol,
                           blocksqp::SymMatrix *&hess, int dmode, int *info) {
    blocksqp::Matrix constrJac;

    constrJac.Dimension(nCon, nVar).Initialize(inf);
    evaluate(xi, lambda, objval, constr, gradObj, constrJac, hess, dmode, info);

    // Convert to sparse format
    if (dmode != 0)
      convertJacobian(constrJac, jacNz, jacIndRow, jacIndCol, 0);
  }


  void BlocksqpProblem::evaluate(const blocksqp::Matrix &xi, const blocksqp::Matrix &lambda,
                           double *objval, blocksqp::Matrix &constr,
                           blocksqp::Matrix &gradObj, blocksqp::Matrix &constrJac,
                           blocksqp::SymMatrix *&hess,
                           int dmode, int *info) {
    *info = 0;

    /*
     * min   x1**2 - 0.5*x2**2
     * s.t.  x1 - x2 = 0
     */
    if (dmode >= 0) {
      *objval = xi(0)*xi(0) - 0.5*xi(1)*xi(1);
      constr(0) = xi(0) - xi(1);
    }

    if (dmode > 0) {
      gradObj(0) = 2.0 * xi(0);
      gradObj(1) = -xi(1);

      constrJac(0, 0) = 1.0;
      constrJac(0, 1) = -1.0;
    }
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
    create_function("nlp_gf_jg", {"x", "p"}, {"f", "g", "grad:f:x", "jac:g:x"});

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

    // Optimal cost
    if (m->f) *m->f = meth->vars->obj;
    // Get primal solution
    casadi_copy(meth->vars->xi.array, nx_, m->x);
    // Get dual solution
    casadi_copy(meth->vars->lambda.array, nx_, m->lam_x);
    casadi_copy(meth->vars->lambda.array + nx_, ng_, m->lam_g);

    // Clean up
    delete prob;
    delete stats;
    delete opts;
    delete meth;
  }

} // namespace casadi
