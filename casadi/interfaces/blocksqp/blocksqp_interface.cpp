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


  class MyProblem : public blocksqp::Problemspec {
  public:
    blocksqp::Matrix xi0;

  public:
    MyProblem(int nVar_,
              int nCon_,
              int nBlocks_,
              int *BlockIdx_,
              const blocksqp::Matrix &bl_,
              const blocksqp::Matrix &bu_,
              const blocksqp::Matrix &xi0_);

    // Set initial values for xi (and possibly lambda) and parts of the
    // Jacobian that correspond to linear constraints (dense version).
    virtual void initialize(blocksqp::Matrix &xi,
                            blocksqp::Matrix &lambda,
                            blocksqp::Matrix &constrJac);

    // Set initial values for xi (and possibly lambda) and parts of the
    // Jacobian that correspond to linear constraints (sparse version).
    virtual void initialize(blocksqp::Matrix &xi,
                            blocksqp::Matrix &lambda,
                            double *&jacNz,
                            int *&jacIndRow,
                            int *&jacIndCol);

    /// Evaluate objective, constraints, and derivatives (dense version).
    virtual void evaluate(const blocksqp::Matrix &xi,
                          const blocksqp::Matrix &lambda,
                          double *objval,
                          blocksqp::Matrix &constr,
                          blocksqp::Matrix &gradObj,
                          blocksqp::Matrix &constrJac,
                          blocksqp::SymMatrix *&hess,
                          int dmode,
                          int *info);

    /// Evaluate objective, constraints, and derivatives (sparse version).
    virtual void evaluate(const blocksqp::Matrix &xi,
                          const blocksqp::Matrix &lambda,
                          double *objval,
                          blocksqp::Matrix &constr,
                          blocksqp::Matrix &gradObj,
                          double *&jacNz,
                          int *&jacIndRow,
                          int *&jacIndCol,
                          blocksqp::SymMatrix *&hess,
                          int dmode,
                          int *info);

    // Generic method to convert dense constraint Jacobian to a sparse matrix
    // in Harwell--Boeing (column compressed) format.
    virtual void convertJacobian(const blocksqp::Matrix &constrJac,
                                 double *&jacNz,
                                 int *&jacIndRow,
                                 int *&jacIndCol,
                                 bool firstCall = 0);
  };


  MyProblem::MyProblem(int nVar_, int nCon_, int nBlocks_, int *blockIdx_,
                       const blocksqp::Matrix &bl_, const blocksqp::Matrix &bu_,
                       const blocksqp::Matrix &xi0_) {
    nVar = nVar_;
    nCon = nCon_;

    nBlocks = nBlocks_;
    blockIdx = new int[nBlocks+1];
    if (nBlocks == 1) {
      blockIdx[0] = 0;
      blockIdx[1] = nVar;
    } else {
      for (int i=0; i<nBlocks+1; i++) blockIdx[i] = blockIdx_[i];
    }

    bl.Dimension(nVar + nCon).Initialize(-inf);
    bu.Dimension(nVar + nCon).Initialize(inf);

    for (int i=0; i<nVar+nCon; i++) {
      bl(i) = bl_(i);
      bu(i) = bu_(i);
    }

    objLo = -inf;
    objUp = inf;

    xi0.Dimension(nVar);
    for (int i=0; i<nVar; i++)
      xi0(i) = xi0_(i);
  }


  void MyProblem::convertJacobian(const blocksqp::Matrix &constrJac, double *&jacNz,
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


  void MyProblem::initialize(blocksqp::Matrix &xi, blocksqp::Matrix &lambda,
                             blocksqp::Matrix &constrJac) {
    // set initial values for xi and lambda
    lambda.Initialize(0.0);
    for (int i=0; i<nVar; i++)
      xi(i) = xi0(i);
  }


  void MyProblem::initialize(blocksqp::Matrix &xi, blocksqp::Matrix &lambda,
                             double *&jacNz, int *&jacIndRow, int *&jacIndCol) {
    blocksqp::Matrix constrDummy, gradObjDummy, constrJac;
    blocksqp::SymMatrix *hessDummy;
    double objvalDummy;
    int info;

    // set initial values for xi and lambda
    lambda.Initialize(0.0);
    for (int i=0; i<nVar; i++)
      xi(i) = xi0(i);

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

  void MyProblem::evaluate(const blocksqp::Matrix &xi, const blocksqp::Matrix &lambda,
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


  void MyProblem::evaluate(const blocksqp::Matrix &xi, const blocksqp::Matrix &lambda,
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
    MyProblem *prob;
    blocksqp::SQPMethod *meth;
    blocksqp::SQPoptions *opts;
    blocksqp::SQPstats *stats;
    char outpath[255];
    strcpy(outpath, "./");

    // Initial value
    blocksqp::Matrix x0;
    x0.Dimension(nx_);
    for (int i=0; i<nx_; ++i) x0(i) = m->x0 ? m->x0[i] : 0;

    // Bounds on variables and constraints
    blocksqp::Matrix bl, bu;
    bl.Dimension(nx_ + ng_);
    bu.Dimension(nx_ + ng_);
    for (int i=0; i<nx_; ++i) {
      bl(i) = m->lbx ? m->lbx[i] : 0;
      bu(i) = m->ubx ? m->ubx[i] : 0;
    }
    for (int i=0; i<ng_; ++i) {
      bl(nx_ + i) = m->lbg ? m->lbg[i] : 0;
      bu(nx_ + i) = m->ubg ? m->ubg[i] : 0;
    }

    // Variable partition for block Hessian
    int nBlocks = nx_;
    vector<int> blockIdx(nBlocks+1);
    for (int i=0; i<nBlocks+1; i++) blockIdx[i] = i;

    // Create problem evaluation object
    prob = new MyProblem( nx_, ng_, nBlocks, get_ptr(blockIdx), bl, bu, x0 );

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
