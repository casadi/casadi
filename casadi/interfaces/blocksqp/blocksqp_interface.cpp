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
     {{"blocksqp",
       {OT_DICT,
        "Options to be passed to BLOCKSQP"}},
      {"contype",
       {OT_INTVECTOR,
        "Type of constraint"}}
     }
  };

  void Blocksqp::init(const Dict& opts) {
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

    /*------------------------*/
    /* Options for SQP solver */
    /*------------------------*/
    m->param = new blocksqp::SQPoptions();
    m->param->opttol = 1.0e-12;
    m->param->nlinfeastol = 1.0e-12;

    // 0: no globalization, 1: filter line search
    m->param->globalization = 0;
    // 0: (scaled) identity, 1: SR1, 2: BFGS
    m->param->hessUpdate = 0;
    // 0: initial Hessian is diagonal matrix, 1: scale initial Hessian according to Nocedal p.143,
    // 2: scale initial Hessian with Oren-Luenberger factor 3: scale initial Hessian
    //    with geometric mean of 1 and 2
    // 4: scale Hessian in every step with centered Oren-Luenberger sizing according to Tapia paper
    m->param->hessScaling = 0;
    // scaling strategy for fallback BFGS update if SR1 and globalization is used
    m->param->fallbackScaling = 0;
    // Size of limited memory
    m->param->hessLimMem = 0;
    // If too many updates are skipped, reset Hessian
    m->param->maxConsecSkippedUpdates = 200;
    // 0: full space Hessian approximation (ignore block structure), 1: blockwise updates
    m->param->blockHess = 0;
    m->param->whichSecondDerv = 0;
    m->param->sparseQP = 1;
    m->param->printLevel = 2;


    /*-------------------------------------------------*/
    /* Create blockSQP method object and run algorithm */
    /*-------------------------------------------------*/
    m->stats = new blocksqp::SQPstats(outpath);

    // Check if there are options that are infeasible and set defaults accordingly
    m->param->optionsConsistency();

    m->vars = new blocksqp::SQPiterate(m->prob, m->param, 1);

    if (m->param->sparseQP < 2) {
      m->qp = new qpOASES::SQProblem(m->prob->nVar, m->prob->nCon);
      m->qpSave = new qpOASES::SQProblem(m->prob->nVar, m->prob->nCon);
    } else {
      m->qp = new qpOASES::SQProblemSchur(m->prob->nVar, m->prob->nCon, qpOASES::HST_UNKNOWN, 50);
      m->qpSave = new qpOASES::SQProblemSchur(m->prob->nVar,
        m->prob->nCon, qpOASES::HST_UNKNOWN, 50);
    }

    m->initCalled = false;

    // Print header and information about the algorithmic parameters
    printInfo(m, m->param->printLevel);

    // Open output files
    m->stats->initStats(m->param);
    m->vars->initIterate(m->param);

    // Initialize filter with pair (maxConstrViolation, objLowerBound)
    initializeFilter(m);

    // Set initial values for all xi and set the Jacobian for linear constraints
    if (m->param->sparseQP) {
      m->prob->initialize(m->vars->xi, m->vars->lambda, m->vars->jacNz,
        m->vars->jacIndRow, m->vars->jacIndCol);
    } else {
      m->prob->initialize(m->vars->xi, m->vars->lambda, m->vars->constrJac);
    }

    m->initCalled = true;




    ret = run(m, 100);
    finish(m);
    if (ret==1) casadi_warning("Maximum number of iterations reached");

    // Get optimal cost
    if (m->f) *m->f = m->vars->obj;
    // Get primal solution
    casadi_copy(m->vars->xi.array, nx_, m->x);
    // Get dual solution (simple bounds)
    if (m->lam_x) {
      casadi_copy(m->vars->lambda.array, nx_, m->lam_x);
      casadi_scal(nx_, -1., m->lam_x);
    }
    // Get dual solution (nonlinear bounds)
    if (m->lam_g) {
      casadi_copy(m->vars->lambda.array + nx_, ng_, m->lam_g);
      casadi_scal(ng_, -1., m->lam_g);
    }

    // Clean up
    delete m->prob;
    delete m->stats;
    delete m->param;
    delete m->qp;
    delete m->qpSave;
    delete m->vars;
  }

  int Blocksqp::run(BlocksqpMemory* m, int maxIt, int warmStart) const {
    int it, infoQP = 0, infoEval = 0;
    bool skipLineSearch = false;
    bool hasConverged = false;
    int whichDerv = m->param->whichSecondDerv;

    if (!m->initCalled) {
      printf("init() must be called before run(). Aborting.\n");
      return -1;
    }

    if (warmStart == 0 || m->stats->itCount == 0) {
      // SQP iteration 0

      /// Set initial Hessian approximation
      calcInitialHessian(m);

      /// Evaluate all functions and gradients for xi_0
      if (m->param->sparseQP) {
        m->prob->evaluate(m->vars->xi, m->vars->lambda, &m->vars->obj,
          m->vars->constr, m->vars->gradObj,
                        m->vars->jacNz, m->vars->jacIndRow, m->vars->jacIndCol,
                        m->vars->hess, 1+whichDerv, &infoEval);
      } else {
        m->prob->evaluate(m->vars->xi, m->vars->lambda, &m->vars->obj,
          m->vars->constr, m->vars->gradObj,
                        m->vars->constrJac, m->vars->hess, 1+whichDerv, &infoEval);
      }
      m->stats->nDerCalls++;

      /// Check if converged
      hasConverged = calcOptTol(m);
      m->stats->printProgress(m->prob, m->vars, m->param, hasConverged);
      if (hasConverged)
        return 0;

      m->stats->itCount++;
    }

    /*
     * SQP Loop: during first iteration, m->stats->itCount = 1
     */
    for (it=0; it<maxIt; it++) {
      /// Solve QP subproblem with qpOASES or QPOPT
      updateStepBounds(m, 0);
      infoQP = solveQP(m, m->vars->deltaXi, m->vars->lambdaQP);

      if (infoQP == 1) {
          // 1.) Maximum number of iterations reached
          printf("***Warning! Maximum number of QP iterations exceeded.***\n");
      } else if (infoQP == 2 || infoQP > 3) {
          // 2.) QP error (e.g., unbounded), solve again with pos.def. diagonal matrix (identity)
          printf("***QP error. Solve again with identity matrix.***\n");
          resetHessian(m);
          infoQP = solveQP(m, m->vars->deltaXi, m->vars->lambdaQP);
          if (infoQP) {
            // If there is still an error, terminate.
            printf("***QP error. Stop.***\n");
            return -1;
          } else {
            m->vars->steptype = 1;
          }
        } else if (infoQP == 3) {
        // 3.) QP infeasible, try to restore feasibility
        bool qpError = true;
        skipLineSearch = true; // don't do line search with restoration step

        // Try to reduce constraint violation by heuristic
        if (m->vars->steptype < 2) {
          printf("***QP infeasible. Trying to reduce constraint violation...");
          qpError = feasibilityRestorationHeuristic(m);
          if (!qpError) {
            m->vars->steptype = 2;
            printf("Success.***\n");
          } else {
            printf("Failed.***\n");
          }
        }

        // Invoke feasibility restoration phase
        //if (qpError && m->vars->steptype < 3 && m->param->restoreFeas)
        if (qpError && m->param->restoreFeas && m->vars->cNorm > 0.01 * m->param->nlinfeastol) {
          printf("***Start feasibility restoration phase.***\n");
          m->vars->steptype = 3;
          qpError = feasibilityRestorationPhase(m);
        }

        // If everything failed, abort.
        if (qpError) {
          printf("***QP error. Stop.***\n");
          return -1;
        }
      }

      /// Determine steplength alpha
      if (m->param->globalization == 0 || (m->param->skipFirstGlobalization
        && m->stats->itCount == 1)) {
        // No globalization strategy, but reduce step if function cannot be evaluated
        if (fullstep(m)) {
          printf("***Constraint or objective could not be evaluated at new point. Stop.***\n");
          return -1;
        }
        m->vars->steptype = 0;
      } else if (m->param->globalization == 1 && !skipLineSearch) {
        // Filter line search based on Waechter et al., 2006 (Ipopt paper)
        if (filterLineSearch(m) || m->vars->reducedStepCount > m->param->maxConsecReducedSteps) {
          // Filter line search did not produce a step. Now there are a few things we can try ...
          bool lsError = true;

          // Heuristic 1: Check if the full step reduces the KKT error by at
          // least kappaF, if so, accept the step.
          lsError = kktErrorReduction(m);
          if (!lsError)
            m->vars->steptype = -1;

          // Heuristic 2: Try to reduce constraint violation by closing
          // continuity gaps to produce an admissable iterate
          if (lsError && m->vars->cNorm > 0.01 * m->param->nlinfeastol && m->vars->steptype < 2) {
            // Don't do this twice in a row!

            printf("***Warning! Steplength too short. Trying to reduce constraint violation...");

            // Integration over whole time interval
            lsError = feasibilityRestorationHeuristic(m);
            if (!lsError) {
                m->vars->steptype = 2;
                printf("Success.***\n");
              } else {
              printf("Failed.***\n");
            }
          }

          // Heuristic 3: Recompute step with a diagonal Hessian
          if (lsError && m->vars->steptype != 1 && m->vars->steptype != 2) {
            // After closing continuity gaps, we already take a step with initial Hessian.
            // If this step is not accepted then this will cause an infinite loop!

            printf("***Warning! Steplength too short. "
                  "Trying to find a new step with identity Hessian.***\n");
            m->vars->steptype = 1;

            resetHessian(m);
            continue;
          }

          // If this does not yield a successful step, start restoration phase
          if (lsError && m->vars->cNorm > 0.01 * m->param->nlinfeastol && m->param->restoreFeas) {
            printf("***Warning! Steplength too short. Start feasibility restoration phase.***\n");
            m->vars->steptype = 3;

            // Solve NLP with minimum norm objective
            lsError = feasibilityRestorationPhase(m);
          }

          // If everything failed, abort.
          if (lsError) {
            printf("***Line search error. Stop.***\n");
            return -1;
          }
        } else {
          m->vars->steptype = 0;
        }
      }

      /// Calculate "old" Lagrange gradient: gamma = dL(xi_k, lambda_k+1)
      calcLagrangeGradient(m, m->vars->gamma, 0);

      /// Evaluate functions and gradients at the new xi
      if (m->param->sparseQP) {
        m->prob->evaluate(m->vars->xi, m->vars->lambda, &m->vars->obj, m->vars->constr,
          m->vars->gradObj, m->vars->jacNz, m->vars->jacIndRow,
          m->vars->jacIndCol, m->vars->hess, 1+whichDerv, &infoEval);
      } else {
        m->prob->evaluate(m->vars->xi, m->vars->lambda, &m->vars->obj, m->vars->constr,
            m->vars->gradObj, m->vars->constrJac, m->vars->hess, 1+whichDerv, &infoEval);
      }
      m->stats->nDerCalls++;

      /// Check if converged
      hasConverged = calcOptTol(m);

      /// Print one line of output for the current iteration
      m->stats->printProgress(m->prob, m->vars, m->param, hasConverged);
      if (hasConverged && m->vars->steptype < 2) {
        m->stats->itCount++;
        if (m->param->debugLevel > 2) {
          //printf("Computing finite differences Hessian at the solution ... \n");
          //calcFiniteDiffHessian();
          //m->stats->printHessian(m->prob->nBlocks, m->vars->hess);
          m->stats->dumpQPCpp(m->prob, m->vars, m->qp, m->param->sparseQP);
        }
        return 0; //Convergence achieved!
      }

      /// Calculate difference of old and new Lagrange gradient:
      // gamma = -gamma + dL(xi_k+1, lambda_k+1)
      calcLagrangeGradient(m, m->vars->gamma, 1);

      /// Revise Hessian approximation
      if (m->param->hessUpdate < 4 && !m->param->hessLimMem) {
        calcHessianUpdate(m, m->param->hessUpdate, m->param->hessScaling);
      } else if (m->param->hessUpdate < 4 && m->param->hessLimMem) {
        calcHessianUpdateLimitedMemory(m, m->param->hessUpdate, m->param->hessScaling);
      } else if (m->param->hessUpdate == 4) {
        calcFiniteDiffHessian(m);
      }

      // If limited memory updates  are used, set pointers deltaXi and
      // gamma to the next column in deltaMat and gammaMat
      updateDeltaGamma(m);

      m->stats->itCount++;
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

    m->stats->finish(m->param);
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
    if (m->param->sparseQP) {
      calcLagrangeGradient(m, m->vars->lambda, m->vars->gradObj, m->vars->jacNz,
        m->vars->jacIndRow, m->vars->jacIndCol, gradLagrange, flag);
    } else {
      calcLagrangeGradient(m, m->vars->lambda, m->vars->gradObj, m->vars->constrJac,
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
    calcLagrangeGradient(m, m->vars->gradLagrange, 0);
    m->vars->gradNorm = lInfVectorNorm(m->vars->gradLagrange);
    m->vars->tol = m->vars->gradNorm /(1.0 + lInfVectorNorm(m->vars->lambda));

    // norm of constraint violation
    m->vars->cNorm  = lInfConstraintNorm(m->vars->xi, m->vars->constr, m->prob->bu, m->prob->bl);
    m->vars->cNormS = m->vars->cNorm /(1.0 + lInfVectorNorm(m->vars->xi));

    if (m->vars->tol <= m->param->opttol && m->vars->cNormS <= m->param->nlinfeastol)
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
    if (m->param->sparseQP == 0)
      strcpy(qpString, "dense, reduced Hessian factorization");
    else if (m->param->sparseQP == 1)
      strcpy(qpString, "sparse, reduced Hessian factorization");
    else if (m->param->sparseQP == 2)
      strcpy(qpString, "sparse, Schur complement approach");

    /* Globalization */
    if (m->param->globalization == 0)
      strcpy(globString, "none (full step)");
    else if (m->param->globalization == 1)
      strcpy(globString, "filter line search");

    /* Hessian approximation */
    if (m->param->blockHess && (m->param->hessUpdate == 1 || m->param->hessUpdate == 2))
      strcpy(hessString1, "block ");
    else
      strcpy(hessString1, "");

    if (m->param->hessLimMem && (m->param->hessUpdate == 1 || m->param->hessUpdate == 2))
      strcat(hessString1, "L-");

    /* Fallback Hessian */
    if (m->param->hessUpdate == 1 || m->param->hessUpdate == 4
      || (m->param->hessUpdate == 2 && !m->param->hessDamp)) {
        strcpy(hessString2, hessString1);

        /* Fallback Hessian update type */
        if (m->param->fallbackUpdate == 0) {
          strcat(hessString2, "Id");
        } else if (m->param->fallbackUpdate == 1) {
          strcat(hessString2, "SR1");
        } else if (m->param->fallbackUpdate == 2) {
          strcat(hessString2, "BFGS");
        } else if (m->param->fallbackUpdate == 4) {
          strcat(hessString2, "Finite differences");
        }

        /* Fallback Hessian scaling */
        if (m->param->fallbackScaling == 1) {
          strcat(hessString2, ", SP");
        } else if (m->param->fallbackScaling == 2) {
          strcat(hessString2, ", OL");
        } else if (m->param->fallbackScaling == 3) {
          strcat(hessString2, ", mean");
        } else if (m->param->fallbackScaling == 4) {
          strcat(hessString2, ", selective sizing");
        }
      } else {
        strcpy(hessString2, "-");
    }

    /* First Hessian update type */
    if (m->param->hessUpdate == 0) {
      strcat(hessString1, "Id");
    } else if (m->param->hessUpdate == 1) {
      strcat(hessString1, "SR1");
    } else if (m->param->hessUpdate == 2) {
      strcat(hessString1, "BFGS");
    } else if (m->param->hessUpdate == 4) {
      strcat(hessString1, "Finite differences");
    }

    /* First Hessian scaling */
    if (m->param->hessScaling == 1) {
      strcat(hessString1, ", SP");
    } else if (m->param->hessScaling == 2) {
      strcat(hessString1, ", OL");
    } else if (m->param->hessScaling == 3) {
      strcat(hessString1, ", mean");
    } else if (m->param->hessScaling == 4) {
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
    acceptStep(m, m->vars->deltaXi, m->vars->lambdaQP, alpha, 0);
  }

  void Blocksqp::
  acceptStep(BlocksqpMemory* m, const blocksqp::Matrix &deltaXi,
    const blocksqp::Matrix &lambdaQP, double alpha, int nSOCS) const {
    int k;
    double lStpNorm;

    // Current alpha
    m->vars->alpha = alpha;
    m->vars->nSOCS = nSOCS;

    // Set new xi by accepting the current trial step
    for (k=0; k<m->vars->xi.M(); k++) {
      m->vars->xi(k) = m->vars->trialXi(k);
      m->vars->deltaXi(k) = alpha * deltaXi(k);
    }

    // Store the infinity norm of the multiplier step
    m->vars->lambdaStepNorm = 0.0;
    for (k=0; k<m->vars->lambda.M(); k++)
      if ((lStpNorm = fabs(alpha*lambdaQP(k) - alpha*m->vars->lambda(k))) > m->vars->lambdaStepNorm)
        m->vars->lambdaStepNorm = lStpNorm;

    // Set new multipliers
    for (k=0; k<m->vars->lambda.M(); k++)
      m->vars->lambda(k) = (1.0 - alpha)*m->vars->lambda(k) + alpha*lambdaQP(k);

    // Count consecutive reduced steps
    if (m->vars->alpha < 1.0)
      m->vars->reducedStepCount++;
    else
      m->vars->reducedStepCount = 0;
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
      if (m->prob->bl(nVar+i) != m->param->inf)
        m->vars->deltaBl(nVar+i) = (*alphaSOC)*m->vars->deltaBl(nVar+i) - m->vars->constr(i);
      else
        m->vars->deltaBl(nVar+i) = m->param->inf;

      if (m->prob->bu(nVar+i) != m->param->inf)
        m->vars->deltaBu(nVar+i) = (*alphaSOC)*m->vars->deltaBu(nVar+i) - m->vars->constr(i);
      else
        m->vars->deltaBu(nVar+i) = m->param->inf;
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
        m->vars->trialXi(i) = m->vars->xi(i) + alpha * m->vars->deltaXi(i);

      // Compute problem functions at trial point
      m->prob->evaluate(m->vars->trialXi, &objTrial, m->vars->constr, &info);
      m->stats->nFunCalls++;
      cNormTrial = lInfConstraintNorm(m->vars->trialXi, m->vars->constr, m->prob->bu, m->prob->bl);
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
    cNorm = lInfConstraintNorm(m->vars->xi, m->vars->constr, m->prob->bu, m->prob->bl);

    // Backtracking line search
    for (k=0; k<m->param->maxLineSearch; k++) {
      // Compute new trial point
      for (i=0; i<nVar; i++)
        m->vars->trialXi(i) = m->vars->xi(i) + alpha * m->vars->deltaXi(i);

      // Compute grad(f)^T * deltaXi
      dfTdeltaXi = 0.0;
      for (i=0; i<nVar; i++)
        dfTdeltaXi += m->vars->gradObj(i) * m->vars->deltaXi(i);

      // Compute objective and at ||constr(trialXi)||_1 at trial point
      m->prob->evaluate(m->vars->trialXi, &objTrial, m->vars->constr, &info);
      m->stats->nFunCalls++;
      cNormTrial = lInfConstraintNorm(m->vars->trialXi, m->vars->constr, m->prob->bu, m->prob->bl);
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
      if (cNorm <= m->param->thetaMin) {
        // Switching condition, part 1: grad(f)^T * deltaXi < 0 ?
        if (dfTdeltaXi < 0)
          // Switching condition, part 2: alpha * (- grad(f)^T * deltaXi)**sF
          // > delta * cNorm**sTheta ?
          if (alpha * pow((-dfTdeltaXi), m->param->sF)
          > m->param->delta * pow(cNorm, m->param->sTheta)) {
            // Switching conditions hold: Require satisfaction of Armijo condition for objective
            if (objTrial > m->vars->obj + m->param->eta*alpha*dfTdeltaXi) {
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
      if (cNormTrial < (1.0 - m->param->gammaTheta) * cNorm
      || objTrial < m->vars->obj - m->param->gammaF * cNorm) {
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
    if (k == m->param->maxLineSearch) return 1;

    // Augment the filter if switching condition or Armijo condition does not hold
    if (dfTdeltaXi >= 0) {
      augmentFilter(m, cNormTrial, objTrial);
    } else if (alpha * pow((-dfTdeltaXi), m->param->sF)
      > m->param->delta * pow(cNorm, m->param->sTheta)) {
      // careful with neg. exponents!
      augmentFilter(m, cNormTrial, objTrial);
    } else if (objTrial <= m->vars->obj + m->param->eta*alpha*dfTdeltaXi) {
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

    // m->vars->constr contains result at first trial point: c(xi+deltaXi)
    // m->vars->constrJac, m->vars->AdeltaXi and m->vars->gradObj are unchanged so far.

    // First SOC step
    deltaXiSOC.Dimension(m->vars->deltaXi.M()).Initialize(0.0);
    lambdaQPSOC.Dimension(m->vars->lambdaQP.M()).Initialize(0.0);

    // Second order correction loop
    cNormOld = cNorm;
    for (k=0; k<m->param->maxSOCiter; k++) {
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
        m->vars->trialXi(i) = m->vars->xi(i) + deltaXiSOC(i);
      }

      // Compute objective and ||constr(trialXiSOC)||_1 at SOC trial point
      m->prob->evaluate(m->vars->trialXi, &objTrialSOC, m->vars->constr, &info);
      m->stats->nFunCalls++;
      cNormTrialSOC = lInfConstraintNorm(m->vars->trialXi, m->vars->constr,
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
      if (cNorm <= m->param->thetaMin && swCond) {
        if (objTrialSOC > m->vars->obj + m->param->eta*dfTdeltaXi) {
          // Armijo condition does not hold for SOC step, next SOC step

          // If constraint violation gets worse by SOC, abort
          if (cNormTrialSOC > m->param->kappaSOC * cNormOld) {
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
      if (cNorm > m->param->thetaMin || !swCond) {
        if (cNormTrialSOC < (1.0 - m->param->gammaTheta) * cNorm
        || objTrialSOC < m->vars->obj - m->param->gammaF * cNorm) {
          // found suitable alpha during SOC, stop
          acceptStep(m, deltaXiSOC, lambdaQPSOC, 1.0, nSOCS);
          return true;
        } else {
          // Trial point is dominated by current point, next SOC step

          // If constraint violation gets worse by SOC, abort
          if (cNormTrialSOC > m->param->kappaSOC * cNormOld) {
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
    if (m->param->restoreFeas == 0) return -1;

    casadi_error("not implemented");
    return 0;
  }


  /**
   * Try to (partly) improve constraint violation by satisfying
   * the (pseudo) continuity constraints, i.e. do a single shooting
   * iteration with the current controls and measurement weights q and w
   */
  int Blocksqp::feasibilityRestorationHeuristic(BlocksqpMemory* m) const {
    m->stats->nRestHeurCalls++;

    int info, k;
    double cNormTrial;

    info = 0;

    // Call problem specific heuristic to reduce constraint violation.
    // For shooting methods that means setting consistent values for
    // shooting nodes by one forward integration.
    for (k=0; k<m->prob->nVar; k++) // input: last successful step
      m->vars->trialXi(k) = m->vars->xi(k);
    m->prob->reduceConstrVio(m->vars->trialXi, &info);
    if (info) {
      // If an error occured in restoration heuristics, abort
      return -1;
    }

    // Compute objective and constraints at the new (hopefully feasible) point
    m->prob->evaluate(m->vars->trialXi, &m->vars->obj, m->vars->constr, &info);
    m->stats->nFunCalls++;
    cNormTrial = lInfConstraintNorm(m->vars->trialXi, m->vars->constr, m->prob->bu, m->prob->bl);
    if (info != 0 || m->vars->obj < m->prob->objLo || m->vars->obj > m->prob->objUp
      || !(m->vars->obj == m->vars->obj) || !(cNormTrial == cNormTrial))
      return -1;

    // Is the new point acceptable for the filter?
    if (pairInFilter(m, cNormTrial, m->vars->obj)) {
      // point is in the taboo region, restoration heuristic not successful!
      return -1;
    }

    // If no error occured in the integration all shooting variables now
    // have the values obtained by a single shooting integration.
    // This is done instead of a Newton-like step in the current SQP iteration

    m->vars->alpha = 1.0;
    m->vars->nSOCS = 0;

    // reset reduced step counter
    m->vars->reducedStepCount = 0;

    // Reset lambda
    m->vars->lambda.Initialize(0.0);
    m->vars->lambdaQP.Initialize(0.0);

    // Compute the "step" taken by closing the continuity conditions
    /// \note deltaXi is reset by resetHessian(), so this doesn't matter
    for (k=0; k<m->prob->nVar; k++) {
      //m->vars->deltaXi(k) = m->vars->trialXi(k) - m->vars->xi(k);
      m->vars->xi(k) = m->vars->trialXi(k);
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
      m->vars->trialXi(i) = m->vars->xi(i) + m->vars->deltaXi(i);

    // Compute objective and ||constr(trialXi)|| at trial point
    trialConstr.Dimension(m->prob->nCon).Initialize(0.0);
    m->prob->evaluate(m->vars->trialXi, &objTrial, trialConstr, &info);
    m->stats->nFunCalls++;
    cNormTrial = lInfConstraintNorm(m->vars->trialXi, trialConstr, m->prob->bu, m->prob->bl);
    if (info != 0 || objTrial < m->prob->objLo || objTrial > m->prob->objUp
      || !(objTrial == objTrial) || !(cNormTrial == cNormTrial)) {
      // evaluation error
      return 1;
    }

    // Compute KKT error of the new point

    // scaled norm of Lagrangian gradient
    trialGradLagrange.Dimension(m->prob->nVar).Initialize(0.0);
    if (m->param->sparseQP) {
      calcLagrangeGradient(m, m->vars->lambdaQP, m->vars->gradObj, m->vars->jacNz,
                            m->vars->jacIndRow, m->vars->jacIndCol, trialGradLagrange, 0);
    } else {
      calcLagrangeGradient(m, m->vars->lambdaQP, m->vars->gradObj, m->vars->constrJac,
                            trialGradLagrange, 0);
    }

    trialGradNorm = lInfVectorNorm(trialGradLagrange);
    trialTol = trialGradNorm /(1.0 + lInfVectorNorm(m->vars->lambdaQP));

    if (fmax(cNormTrial, trialTol) < m->param->kappaF * fmax(m->vars->cNorm, m->vars->tol)) {
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
    filter = m->vars->filter;

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
      if ((cNorm >= (1.0 - m->param->gammaTheta) * iter->first ||
           (cNorm < 0.01 * m->param->nlinfeastol && iter->first < 0.01 * m->param->nlinfeastol)) &&
          obj >= iter->second - m->param->gammaF * iter->first) {
        return 1;
      }

    return 0;
  }


  void Blocksqp::initializeFilter(BlocksqpMemory* m) const {
    std::set< std::pair<double, double> >::iterator iter;
    std::pair<double, double> initPair(m->param->thetaMax, m->prob->objLo);

    // Remove all elements
    iter=m->vars->filter->begin();
    while (iter != m->vars->filter->end()) {
      std::set< std::pair<double, double> >::iterator iterToRemove = iter;
      iter++;
      m->vars->filter->erase(iterToRemove);
    }

    // Initialize with pair (maxConstrViolation, objLowerBound);
    m->vars->filter->insert(initPair);
  }


  /**
   * Augment the filter:
   * F_k+1 = F_k U { (c,f) | c > (1-gammaTheta)cNorm and f > obj-gammaF*c
   */
  void Blocksqp::
  augmentFilter(BlocksqpMemory* m, double cNorm, double obj) const {
    std::set< std::pair<double, double> >::iterator iter;
    std::pair<double, double> entry((1.0 - m->param->gammaTheta)*cNorm, obj
      - m->param->gammaF*cNorm);

    // Augment filter by current element
    m->vars->filter->insert(entry);

    // Remove dominated elements
    iter=m->vars->filter->begin();
    while (iter != m->vars->filter->end()) {
      if (iter->first > entry.first && iter->second > entry.second) {
        std::set< std::pair<double, double> >::iterator iterToRemove = iter;
        iter++;
        m->vars->filter->erase(iterToRemove);
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

    for (iBlock=0; iBlock<m->vars->nBlocks; iBlock++)
      //if objective derv is computed exactly, don't set the last block!
      if (!(m->param->whichSecondDerv == 1 && m->param->blockHess
        && iBlock == m->vars->nBlocks-1))
        calcInitialHessian(m, iBlock);
  }


  /**
   * Initial Hessian for one block: Identity matrix
   */
  void Blocksqp::calcInitialHessian(BlocksqpMemory* m, int iBlock) const {
    m->vars->hess[iBlock].Initialize(0.0);

    // Each block is a diagonal matrix
    for (int i=0; i<m->vars->hess[iBlock].M(); i++)
      m->vars->hess[iBlock](i, i) = m->param->iniHessDiag;

    // If we maintain 2 Hessians, also reset the second one
    if (m->vars->hess2 != 0) {
      m->vars->hess2[iBlock].Initialize(0.0);
      for (int i=0; i<m->vars->hess2[iBlock].M(); i++)
        m->vars->hess2[iBlock](i, i) = m->param->iniHessDiag;
    }
  }


  void Blocksqp::resetHessian(BlocksqpMemory* m) const {
    for (int iBlock=0; iBlock<m->vars->nBlocks; iBlock++) {
      //if objective derv is computed exactly, don't set the last block!
      if (!(m->param->whichSecondDerv == 1 && m->param->blockHess
        && iBlock == m->vars->nBlocks - 1))
        resetHessian(m, iBlock);
      }
  }


  void Blocksqp::resetHessian(BlocksqpMemory* m, int iBlock) const {
    blocksqp::Matrix smallDelta, smallGamma;
    int nVarLocal = m->vars->hess[iBlock].M();

    // smallGamma and smallDelta are either subvectors of gamma and delta
    // or submatrices of gammaMat, deltaMat, i.e. subvectors of gamma and delta
    // from m prev. iterations (for L-BFGS)
    smallGamma.Submatrix(m->vars->gammaMat, nVarLocal, m->vars->gammaMat.N(),
      m->vars->blockIdx[iBlock], 0);
    smallDelta.Submatrix(m->vars->deltaMat, nVarLocal, m->vars->deltaMat.N(),
      m->vars->blockIdx[iBlock], 0);

    // Remove past information on Lagrangian gradient difference
    smallGamma.Initialize(0.0);

    // Remove past information on steps
    smallDelta.Initialize(0.0);

    // Remove information on old scalars (used for COL sizing)
    m->vars->deltaNorm(iBlock) = 1.0;
    m->vars->deltaGamma(iBlock) = 0.0;
    m->vars->deltaNormOld(iBlock) = 1.0;
    m->vars->deltaGammaOld(iBlock) = 0.0;

    m->vars->noUpdateCounter[iBlock] = -1;

    calcInitialHessian(m, iBlock);
  }

  /**
   * Approximate Hessian by finite differences
   */
  int Blocksqp::calcFiniteDiffHessian(BlocksqpMemory* m) const {
    int iVar, jVar, k, iBlock, maxBlock, info, idx, idx1, idx2;
    double dummy, lowerVio, upperVio;
    blocksqp::Matrix pert;
    blocksqp::SQPiterate varsP = blocksqp::SQPiterate(*m->vars);

    const double myDelta = 1.0e-4;
    const double minDelta = 1.0e-6;

    pert.Dimension(m->prob->nVar);

    info = 0;

    // Find out the largest block
    maxBlock = 0;
    for (iBlock=0; iBlock<m->vars->nBlocks; iBlock++)
      if (m->vars->blockIdx[iBlock+1] - m->vars->blockIdx[iBlock] > maxBlock)
        maxBlock = m->vars->blockIdx[iBlock+1] - m->vars->blockIdx[iBlock];

    // Compute original Lagrange gradient
    calcLagrangeGradient(m, m->vars->lambda, m->vars->gradObj, m->vars->jacNz,
      m->vars->jacIndRow, m->vars->jacIndCol, m->vars->gradLagrange, 0);

    for (iVar = 0; iVar<maxBlock; iVar++) {
      pert.Initialize(0.0);

      // Perturb all blocks simultaneously
      for (iBlock=0; iBlock<m->vars->nBlocks; iBlock++) {
        idx = m->vars->blockIdx[iBlock] + iVar;
        // Skip blocks that have less than iVar variables
        if (idx < m->vars->blockIdx[iBlock+1]) {
          pert(idx) = myDelta * fabs(m->vars->xi(idx));
          pert(idx) = fmax(pert(idx), minDelta);

          // If perturbation violates upper bound, try to perturb with negative
          upperVio = m->vars->xi(idx) + pert(idx) - m->prob->bu(idx);
          if (upperVio > 0) {
            lowerVio = m->prob->bl(idx) -  (m->vars->xi(idx) - pert(idx));
            // If perturbation violates also lower bound, take the largest perturbation possible
            if (lowerVio > 0) {
              if (lowerVio > upperVio)
                pert(idx) = -lowerVio;
              else
                pert(idx) = upperVio;
            } else {
              // If perturbation does not violate lower bound, take -computed perturbation
              pert(idx) = -pert(idx);
            }
          }
        }
      }

      // Add perturbation
      for (k=0; k<m->prob->nVar; k++)
        m->vars->xi(k) += pert(k);

      // Compute perturbed Lagrange gradient
      if (m->param->sparseQP) {
        m->prob->evaluate(m->vars->xi, m->vars->lambda, &dummy, varsP.constr, varsP.gradObj,
                        varsP.jacNz, varsP.jacIndRow, varsP.jacIndCol, m->vars->hess, 1, &info);
        calcLagrangeGradient(m, m->vars->lambda, varsP.gradObj, varsP.jacNz,
          varsP.jacIndRow, varsP.jacIndCol, varsP.gradLagrange, 0);
      } else {
        m->prob->evaluate(m->vars->xi, m->vars->lambda, &dummy, varsP.constr, varsP.gradObj,
          varsP.constrJac, m->vars->hess, 1, &info);
        calcLagrangeGradient(m, m->vars->lambda, varsP.gradObj,
          varsP.constrJac, varsP.gradLagrange, 0);
      }

      // Compute finite difference approximations: one column in every block
      for (iBlock=0; iBlock<m->vars->nBlocks; iBlock++) {
        idx1 = m->vars->blockIdx[iBlock] + iVar;
        // Skip blocks that have less than iVar variables
        if (idx1 < m->vars->blockIdx[iBlock+1]) {
          for (jVar=iVar; jVar<m->vars->blockIdx[iBlock+1]-m->vars->blockIdx[iBlock]; jVar++) {
            // Take symmetrized matrices
              idx2 = m->vars->blockIdx[iBlock] + jVar;
              m->vars->hess[iBlock](iVar, jVar) = varsP.gradLagrange(idx1)
                - m->vars->gradLagrange(idx2);
              m->vars->hess[iBlock](iVar, jVar) += (varsP.gradLagrange(idx2)
                - m->vars->gradLagrange(idx1));
              m->vars->hess[iBlock](iVar, jVar) *= 0.5 / pert(idx1);
            }
        }
      }

      // Subtract perturbation
      for (k=0; k<m->prob->nVar; k++) m->vars->xi(k) -= pert(k);
    }

    return info;
  }


  void Blocksqp::sizeInitialHessian(BlocksqpMemory* m, const blocksqp::Matrix &gamma,
    const blocksqp::Matrix &delta, int iBlock, int option) const {
    int i, j;
    double scale;
    double myEps = 1.0e3 * m->param->eps;

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
      for (i=0; i<m->vars->hess[iBlock].M(); i++)
        for (j=i; j<m->vars->hess[iBlock].M(); j++)
          m->vars->hess[iBlock](i, j) *= scale;
    } else {
      scale = 1.0;
    }

    // statistics: average sizing factor
    m->stats->averageSizingFactor += scale;
  }


  void Blocksqp::sizeHessianCOL(BlocksqpMemory* m, const blocksqp::Matrix &gamma,
    const blocksqp::Matrix &delta, int iBlock) const {
    int i, j;
    double theta, scale, myEps = 1.0e3 * m->param->eps;
    double deltaNorm, deltaNormOld, deltaGamma, deltaGammaOld, deltaBdelta;

    // Get sTs, sTs_, sTy, sTy_, sTBs
    deltaNorm = m->vars->deltaNorm(iBlock);
    deltaGamma = m->vars->deltaGamma(iBlock);
    deltaNormOld = m->vars->deltaNormOld(iBlock);
    deltaGammaOld = m->vars->deltaGammaOld(iBlock);
    deltaBdelta = 0.0;
    for (i=0; i<delta.M(); i++)
      for (j=0; j<delta.M(); j++)
        deltaBdelta += delta(i) * m->vars->hess[iBlock](i, j) * delta(j);

    // Centered Oren-Luenberger factor
    if (m->vars->noUpdateCounter[iBlock] == -1) {
      // in the first iteration, this should equal the OL factor
      theta = 1.0;
    } else {
      theta = fmin(m->param->colTau1, m->param->colTau2 * deltaNorm);
    }
    if (deltaNorm > myEps && deltaNormOld > myEps) {
      scale = (1.0 - theta)*deltaGammaOld / deltaNormOld + theta*deltaBdelta / deltaNorm;
      if (scale > m->param->eps)
        scale = ((1.0 - theta)*deltaGammaOld / deltaNormOld + theta*deltaGamma / deltaNorm) / scale;
    } else {
      scale = 1.0;
    }

    // Size only if factor is between zero and one
    if (scale < 1.0 && scale > 0.0) {
      scale = fmax(m->param->colEps, scale);
      //printf("Sizing value (COL) block %i = %g\n", iBlock, scale);
      for (i=0; i<m->vars->hess[iBlock].M(); i++)
        for (j=i; j<m->vars->hess[iBlock].M(); j++)
          m->vars->hess[iBlock](i, j) *= scale;

      // statistics: average sizing factor
      m->stats->averageSizingFactor += scale;
    } else {
      m->stats->averageSizingFactor += 1.0;
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
    if (m->param->whichSecondDerv == 1 && m->param->blockHess)
      nBlocks = m->vars->nBlocks - 1;
    else
      nBlocks = m->vars->nBlocks;

    // Statistics: how often is damping active, what is the average COL sizing factor?
    m->stats->hessDamped = 0;
    m->stats->averageSizingFactor = 0.0;

    for (iBlock=0; iBlock<nBlocks; iBlock++) {
      nVarLocal = m->vars->hess[iBlock].M();

      // smallGamma and smallDelta are subvectors of gamma and delta,
      // corresponding to partially separability
      smallGamma.Submatrix(m->vars->gammaMat, nVarLocal, m->vars->gammaMat.N(),
        m->vars->blockIdx[iBlock], 0);
      smallDelta.Submatrix(m->vars->deltaMat, nVarLocal, m->vars->deltaMat.N(),
        m->vars->blockIdx[iBlock], 0);

      // Is this the first iteration or the first after a Hessian reset?
      firstIter = (m->vars->noUpdateCounter[iBlock] == -1);

      // Update sTs, sTs_ and sTy, sTy_
      m->vars->deltaNormOld(iBlock) = m->vars->deltaNorm(iBlock);
      m->vars->deltaGammaOld(iBlock) = m->vars->deltaGamma(iBlock);
      m->vars->deltaNorm(iBlock) = adotb(smallDelta, smallDelta);
      m->vars->deltaGamma(iBlock) = adotb(smallDelta, smallGamma);

      // Sizing before the update
      if (hessScaling < 4 && firstIter)
        sizeInitialHessian(m, smallGamma, smallDelta, iBlock, hessScaling);
      else if (hessScaling == 4)
        sizeHessianCOL(m, smallGamma, smallDelta, iBlock);

      // Compute the new update
      if (updateType == 1) {
        calcSR1(m, smallGamma, smallDelta, iBlock);

        // Prepare to compute fallback update as well
        m->vars->hess = m->vars->hess2;

        // Sizing the fallback update
        if (m->param->fallbackScaling < 4 && firstIter)
          sizeInitialHessian(m, smallGamma, smallDelta, iBlock, m->param->fallbackScaling);
        else if (m->param->fallbackScaling == 4)
          sizeHessianCOL(m, smallGamma, smallDelta, iBlock);

        // Compute fallback update
        if (m->param->fallbackUpdate == 2)
          calcBFGS(m, smallGamma, smallDelta, iBlock);

        // Reset pointer
        m->vars->hess = m->vars->hess1;
      } else if (updateType == 2) {
        calcBFGS(m, smallGamma, smallDelta, iBlock);
      }

      // If an update is skipped to often, reset Hessian block
      if (m->vars->noUpdateCounter[iBlock] > m->param->maxConsecSkippedUpdates) {
        resetHessian(m, iBlock);
      }
    }

    // statistics: average sizing factor
    m->stats->averageSizingFactor /= nBlocks;
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
    if (m->param->whichSecondDerv == 1 && m->param->blockHess) {
      nBlocks = m->vars->nBlocks - 1;
    } else {
      nBlocks = m->vars->nBlocks;
    }

    // Statistics: how often is damping active, what is the average COL sizing factor?
    m->stats->hessDamped = 0;
    m->stats->hessSkipped = 0;
    m->stats->averageSizingFactor = 0.0;

    for (iBlock=0; iBlock<nBlocks; iBlock++) {
      nVarLocal = m->vars->hess[iBlock].M();

      // smallGamma and smallDelta are submatrices of gammaMat, deltaMat,
      // i.e. subvectors of gamma and delta from m prev. iterations
      smallGamma.Submatrix(m->vars->gammaMat, nVarLocal, m->vars->gammaMat.N(),
        m->vars->blockIdx[iBlock], 0);
      smallDelta.Submatrix(m->vars->deltaMat, nVarLocal, m->vars->deltaMat.N(),
        m->vars->blockIdx[iBlock], 0);

      // Memory structure
      if (m->stats->itCount > smallGamma.N()) {
        m2 = smallGamma.N();
        posOldest = m->stats->itCount % m2;
        posNewest = (m->stats->itCount-1) % m2;
      } else {
        m2 = m->stats->itCount;
        posOldest = 0;
        posNewest = m2-1;
      }

      // Set B_0 (pretend it's the first step)
      calcInitialHessian(m, iBlock);
      m->vars->deltaNorm(iBlock) = 1.0;
      m->vars->deltaNormOld(iBlock) = 1.0;
      m->vars->deltaGamma(iBlock) = 0.0;
      m->vars->deltaGammaOld(iBlock) = 0.0;
      m->vars->noUpdateCounter[iBlock] = -1;

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
        m->vars->deltaNormOld(iBlock) = m->vars->deltaNorm(iBlock);
        m->vars->deltaGammaOld(iBlock) = m->vars->deltaGamma(iBlock);
        m->vars->deltaNorm(iBlock) = adotb(deltai, deltai);
        m->vars->deltaGamma(iBlock) = adotb(gammai, deltai);

        // Save statistics, we want to record them only for the most recent update
        averageSizingFactor = m->stats->averageSizingFactor;
        hessDamped = m->stats->hessDamped;
        hessSkipped = m->stats->hessSkipped;

        // Selective sizing before the update
        if (hessScaling == 4) sizeHessianCOL(m, gammai, deltai, iBlock);

        // Compute the new update
        if (updateType == 1) {
          calcSR1(m, gammai, deltai, iBlock);
        } else if (updateType == 2) {
          calcBFGS(m, gammai, deltai, iBlock);
        }

        m->stats->nTotalUpdates++;

        // Count damping statistics only for the most recent update
        if (pos != posNewest) {
          m->stats->hessDamped = hessDamped;
          m->stats->hessSkipped = hessSkipped;
          if (hessScaling == 4)
            m->stats->averageSizingFactor = averageSizingFactor;
        }
      }

      // If an update is skipped to often, reset Hessian block
      if (m->vars->noUpdateCounter[iBlock] > m->param->maxConsecSkippedUpdates) {
        resetHessian(m, iBlock);
      }
    }
    //blocks
    m->stats->averageSizingFactor /= nBlocks;
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
     * Note that m->vars->gamma needs to remain unchanged!
     * This may be important in a limited memory context:
     * When information is "forgotten", B_i-1 is different and the
     *  original gamma might lead to an undamped update with the new B_i-1! */
    blocksqp::Matrix gamma2 = gamma;

    B = &m->vars->hess[iBlock];

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
    h2 = m->vars->deltaGamma(iBlock);

    /* Powell's damping strategy to maintain pos. def. (Nocedal/Wright p.537; SNOPT paper)
     * Interpolates between current approximation and unmodified BFGS */
    damped = 0;
    if (m->param->hessDamp)
      if (h2 < m->param->hessDampFac * h1 / m->vars->alpha && fabs(h1 - h2) > 1.0e-12) {
        // At the first iteration h1 and h2 are equal due to COL scaling

        thetaPowell = (1.0 - m->param->hessDampFac)*h1 / (h1 - h2);

        // Redefine gamma and h2 = delta^T * gamma
        h2 = 0.0;
        for (i=0; i<dim; i++) {
          gamma2(i) = thetaPowell*gamma2(i) + (1.0 - thetaPowell)*Bdelta(i);
          h2 += delta(i) * gamma2(i);
        }

        // Also redefine deltaGamma for computation of sizing factor in the next iteration
        m->vars->deltaGamma(iBlock) = h2;

        damped = 1;
      }

    // For statistics: count number of damped blocks
    m->stats->hessDamped += damped;

    // B_k+1 = B_k - Bdelta * (Bdelta)^T / h1 + gamma * gamma^T / h2
    double myEps = 1.0e2 * m->param->eps;
    if (fabs(h1) < myEps || fabs(h2) < myEps) {
      // don't perform update because of bad condition, might introduce negative eigenvalues
      m->vars->noUpdateCounter[iBlock]++;
      m->stats->hessDamped -= damped;
      m->stats->hessSkipped++;
      m->stats->nTotalSkippedUpdates++;
    } else {
      for (i=0; i<dim; i++)
        for (j=i; j<dim; j++)
          (*B)(i, j) = (*B)(i, j) - Bdelta(i) * Bdelta(j) / h1
            + gamma2(i) * gamma2(j) / h2;

      m->vars->noUpdateCounter[iBlock] = 0;
    }
  }


  void Blocksqp::
  calcSR1(BlocksqpMemory* m, const blocksqp::Matrix &gamma, const blocksqp::Matrix &delta,
    int iBlock) const {
    int i, j, k, dim = gamma.M();
    blocksqp::Matrix gmBdelta;
    blocksqp::SymMatrix *B;
    double myEps = 1.0e2 * m->param->eps;
    double r = 1.0e-8;
    double h = 0.0;

    B = &m->vars->hess[iBlock];

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
      m->vars->noUpdateCounter[iBlock]++;
      m->stats->hessSkipped++;
      m->stats->nTotalSkippedUpdates++;
    } else {
      for (i=0; i<dim; i++)
        for (j=i; j<dim; j++)
          (*B)(i, j) = (*B)(i, j) + gmBdelta(i) * gmBdelta(j) / h;
      m->vars->noUpdateCounter[iBlock] = 0;
    }
  }


  /**
   * Set deltaXi and gamma as a column in the matrix containing
   * the m most recent delta and gamma
   */
  void Blocksqp::updateDeltaGamma(BlocksqpMemory* m) const {
    int nVar = m->vars->gammaMat.M();
    int m2 = m->vars->gammaMat.N();

    if (m2 == 1)
      return;

    m->vars->deltaXi.Submatrix(m->vars->deltaMat, nVar, 1, 0, m->stats->itCount % m2);
    m->vars->gamma.Submatrix(m->vars->gammaMat, nVar, 1, 0, m->stats->itCount % m2);
  }

  void Blocksqp::
  computeNextHessian(BlocksqpMemory* m, int idx, int maxQP) const {
    // Compute fallback update only once
    if (idx == 1) {
        // Switch storage
        m->vars->hess = m->vars->hess2;

        // If last block contains exact Hessian, we need to copy it
        if (m->param->whichSecondDerv == 1)
          for (int i=0; i<m->vars->hess[m->prob->nBlocks-1].M(); i++)
            for (int j=i; j<m->vars->hess[m->prob->nBlocks-1].N(); j++)
              m->vars->hess2[m->prob->nBlocks-1](i, j) = m->vars->hess1[m->prob->nBlocks-1](i, j);

        // Limited memory: compute fallback update only when needed
        if (m->param->hessLimMem) {
            m->stats->itCount--;
            int hessDampSave = m->param->hessDamp;
            m->param->hessDamp = 1;
            calcHessianUpdateLimitedMemory(m, m->param->fallbackUpdate, m->param->fallbackScaling);
            m->param->hessDamp = hessDampSave;
            m->stats->itCount++;
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
        for (int iBlock=0; iBlock<m->vars->nBlocks; iBlock++)
          for (int i=0; i<m->vars->hess[iBlock].M(); i++)
            for (int j=i; j<m->vars->hess[iBlock].N(); j++) {
                m->vars->hess2[iBlock](i, j) *= mu;
                m->vars->hess2[iBlock](i, j) += mu1 * m->vars->hess1[iBlock](i, j);
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
    if (m->param->globalization == 1 &&
        m->param->hessUpdate == 1 &&
        matricesChanged &&
        m->stats->itCount > 1) {
        maxQP = m->param->maxConvQP + 1;
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
        if (m->param->sparseQP) {
            A = new qpOASES::SparseMatrix(m->prob->nCon, m->prob->nVar,
                                           m->vars->jacIndRow, m->vars->jacIndCol, m->vars->jacNz);
          } else {
            // transpose Jacobian (qpOASES needs row major arrays)
            Transpose(m->vars->constrJac, jacT);
            A = new qpOASES::DenseMatrix(m->prob->nCon, m->prob->nVar, m->prob->nVar, jacT.ARRAY());
          }
      }
    double *g = m->vars->gradObj.ARRAY();
    double *lb = m->vars->deltaBl.ARRAY();
    double *lu = m->vars->deltaBu.ARRAY();
    double *lbA = m->vars->deltaBl.ARRAY() + m->prob->nVar;
    double *luA = m->vars->deltaBu.ARRAY() + m->prob->nVar;

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
    double cpuTime = matricesChanged ? m->param->maxTimeQP : 0.1*m->param->maxTimeQP;
    int maxIt = matricesChanged ? m->param->maxItQP : 0.1*m->param->maxItQP;
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
          m->stats->qpResolve++;
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
            if (m->param->sparseQP) {
                // Convert block-Hessian to sparse format
                m->vars->convertHessian(m->prob, m->param->eps, m->vars->hess, m->vars->hessNz,
                                      m->vars->hessIndRow, m->vars->hessIndCol, m->vars->hessIndLo);
                H = new qpOASES::SymSparseMat(m->prob->nVar, m->prob->nVar,
                                               m->vars->hessIndRow, m->vars->hessIndCol,
                                               m->vars->hessNz);
                dynamic_cast<qpOASES::SymSparseMat*>(H)->createDiagInfo();
              } else {
                // Convert block-Hessian to double array
                m->vars->convertHessian(m->prob, m->param->eps, m->vars->hess);
                H = new qpOASES::SymDenseMat(m->prob->nVar, m->prob->nVar,
                  m->prob->nVar, m->vars->hessNz);
              }
          }

        /*
         * Call qpOASES
         */
        if (m->param->debugLevel > 2) {
          m->stats->dumpQPCpp(m->prob, m->vars, m->qp, m->param->sparseQP);
        }
        if (matricesChanged) {
            maxIt = m->param->maxItQP;
            cpuTime = m->param->maxTimeQP;
            if (m->qp->getStatus() == qpOASES::QPS_HOMOTOPYQPSOLVED ||
                m->qp->getStatus() == qpOASES::QPS_SOLVED) {
                ret = m->qp->hotstart(H, g, A, lb, lu, lbA, luA, maxIt, &cpuTime);
              } else {
                ret = m->qp->init(H, g, A, lb, lu, lbA, luA, maxIt, &cpuTime);
              }
          } else if (!matricesChanged) {
            // Second order correction: H and A do not change
            maxIt = 0.1*m->param->maxItQP;
            cpuTime = 0.1*m->param->maxTimeQP;
            ret = m->qp->hotstart(g, lb, lu, lbA, luA, maxIt, &cpuTime);
          }

        /*
         * Check assumption (G3*) if nonconvex QP was solved
         */
        if (l < maxQP-1 && matricesChanged) {
            if (ret == qpOASES::SUCCESSFUL_RETURN) {
                if (m->param->sparseQP == 2) {
                  ret = solAna.checkCurvatureOnStronglyActiveConstraints(
                    dynamic_cast<qpOASES::SQProblemSchur*>(m->qp));
                } else {
                  ret = solAna.checkCurvatureOnStronglyActiveConstraints(m->qp);
                }
              }

            if (ret == qpOASES::SUCCESSFUL_RETURN) {
              // QP was solved successfully and curvature is positive after removing bounds
                m->stats->qpIterations = maxIt + 1;
                break; // Success!
              } else {
              // QP solution is rejected, save statistics
                if (ret == qpOASES::RET_SETUP_AUXILIARYQP_FAILED)
                  m->stats->qpIterations2++;
                else
                  m->stats->qpIterations2 += maxIt + 1;
                m->stats->rejectedSR1++;
              }
          } else {
            // Convex QP was solved, no need to check assumption (G3*)
            m->stats->qpIterations += maxIt + 1;
          }

      } // End of QP solving loop

    /*
     * Post-processing
     */

    // Get solution from qpOASES
    m->qp->getPrimalSolution(deltaXi.ARRAY());
    m->qp->getDualSolution(lambdaQP.ARRAY());
    m->vars->qpObj = m->qp->getObjVal();

    // Compute constrJac*deltaXi, need this for second order correction step
    if (m->param->sparseQP) {
      Atimesb(m->vars->jacNz, m->vars->jacIndRow, m->vars->jacIndCol, deltaXi, m->vars->AdeltaXi);
    } else {
      Atimesb(m->vars->constrJac, deltaXi, m->vars->AdeltaXi);
    }

    // Print qpOASES error code, if any
    if (ret != qpOASES::SUCCESSFUL_RETURN && matricesChanged)
      printf("qpOASES error message: \"%s\"\n",
              qpOASES::getGlobalMessageHandler()->getErrorCodeMessage(ret));

    // Point Hessian again to the first Hessian
    m->vars->hess = m->vars->hess1;

    /* For full-memory Hessian: Restore fallback Hessian if convex combinations
     * were used during the loop */
    if (!m->param->hessLimMem && maxQP > 2 && matricesChanged) {
        double mu = 1.0 / l;
        double mu1 = 1.0 - mu;
        int nBlocks = (m->param->whichSecondDerv == 1) ? m->vars->nBlocks-1 : m->vars->nBlocks;
        for (int iBlock=0; iBlock<nBlocks; iBlock++)
          for (int i=0; i<m->vars->hess[iBlock].M(); i++)
            for (int j=i; j<m->vars->hess[iBlock].N(); j++) {
                m->vars->hess2[iBlock](i, j) *= mu;
                m->vars->hess2[iBlock](i, j) += mu1 * m->vars->hess1[iBlock](i, j);
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
      if (m->prob->bl(i) != m->param->inf)
        m->vars->deltaBl(i) = m->prob->bl(i) - m->vars->xi(i);
      else
        m->vars->deltaBl(i) = m->param->inf;

      if (m->prob->bu(i) != m->param->inf)
        m->vars->deltaBu(i) = m->prob->bu(i) - m->vars->xi(i);
      else
        m->vars->deltaBu(i) = m->param->inf;
    }

    // Bounds on linearized constraints
    for (i=0; i<nCon; i++) {
      if (m->prob->bl(nVar+i) != m->param->inf) {
        m->vars->deltaBl(nVar+i) = m->prob->bl(nVar+i) - m->vars->constr(i);
        if (soc) m->vars->deltaBl(nVar+i) += m->vars->AdeltaXi(i);
      } else {
        m->vars->deltaBl(nVar+i) = m->param->inf;
      }

      if (m->prob->bu(nVar+i) != m->param->inf) {
        m->vars->deltaBu(nVar+i) = m->prob->bu(nVar+i) - m->vars->constr(i);
        if (soc) m->vars->deltaBu(nVar+i) += m->vars->AdeltaXi(i);
      } else {
        m->vars->deltaBu(nVar+i) = m->param->inf;
      }
    }
  }

} // namespace casadi
