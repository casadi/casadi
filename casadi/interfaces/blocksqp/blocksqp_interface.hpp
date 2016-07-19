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


#ifndef CASADI_BLOCKSQP_INTERFACE_HPP
#define CASADI_BLOCKSQP_INTERFACE_HPP

#include <casadi/interfaces/blocksqp/casadi_nlpsol_blocksqp_export.h>
#include <blocksqp.hpp>
#include "casadi/core/function/nlpsol_impl.hpp"

/** \defgroup plugin_Nlpsol_blocksqp
  * This is a modified version of blockSQP by Janka et al.
  *
  * \author Dennis Janka, Joel Andersson
  * \date 2012-2015, 2016
*/

/** \pluginsection{Nlpsol,blocksqp} */

/// \cond INTERNAL
namespace casadi {
  // Forward declaration
  class BlocksqpInterface;

  struct CASADI_NLPSOL_BLOCKSQP_EXPORT BlocksqpMemory : public NlpsolMemory {
    blocksqp::Problemspec* prob;
    blocksqp::SQPiterate* vars;
    blocksqp::SQPoptions* param;
    blocksqp::SQPstats* stats;
    qpOASES::SQProblem* qp;
    qpOASES::SQProblem* qpSave;
    bool initCalled;

    double* jac;
  };

  // Problem class
  class BlocksqpProblem : public blocksqp::Problemspec {
  public:
    const BlocksqpInterface& self;
    BlocksqpMemory* m;
    BlocksqpProblem(const BlocksqpInterface& self, BlocksqpMemory* m);

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
  };


  /** \brief \pluginbrief{Nlpsol,blocksqp}
     @copydoc Nlpsol_doc
     @copydoc plugin_Nlpsol_blocksqp
  */
  class CASADI_NLPSOL_BLOCKSQP_EXPORT BlocksqpInterface : public Nlpsol {
  public:
    explicit BlocksqpInterface(const std::string& name, const Function& nlp);
    virtual ~BlocksqpInterface();

    // Get name of the plugin
    virtual const char* plugin_name() const { return "blocksqp";}

    /** \brief  Create a new NLP Solver */
    static Nlpsol* creator(const std::string& name, const Function& nlp) {
      return new BlocksqpInterface(name, nlp);
    }

    ///@{
    /** \brief Options */
    static Options options_;
    virtual const Options& get_options() const { return options_;}
    ///@}

    // Initialize the solver
    virtual void init(const Dict& opts);

    /** \brief Create memory block */
    virtual void* alloc_memory() const { return new BlocksqpMemory();}

    /** \brief Free memory block */
    virtual void free_memory(void *mem) const { delete static_cast<BlocksqpMemory*>(mem);}

    /** \brief Initalize memory block */
    virtual void init_memory(void* mem) const;

    /** \brief Set the (persistent) work vectors */
    virtual void set_work(void* mem, const double**& arg, double**& res,
                          int*& iw, double*& w) const;

    // Solve the NLP
    virtual void solve(void* mem) const;

    /// A documentation string
    static const std::string meta_doc;

    // Block partitioning
    std::vector<int> blocks_;

    // Jacobian sparsity
    Sparsity sp_jac_;

    // TO BE REFACTORED

    /// Main Loop of SQP method
    int run(BlocksqpMemory* m, int maxIt, int warmStart = 0) const;
    /// Call after the last call of run, to close output files etc.
    void finish(BlocksqpMemory* m) const;
    /// Print information about the SQP method
    void printInfo(BlocksqpMemory* m, int printLevel) const;
    /// Compute gradient of Lagrangian function (dense version)
    void calcLagrangeGradient(BlocksqpMemory* m, const blocksqp::Matrix &lambda,
      const blocksqp::Matrix &gradObj,
      const blocksqp::Matrix &constrJacFull, blocksqp::Matrix &gradLagrange, int flag) const;
    /// Compute gradient of Lagrangian function (sparse version)
    void calcLagrangeGradient(BlocksqpMemory* m, const blocksqp::Matrix &lambda,
      const blocksqp::Matrix &gradObj, double *jacNz, int *jacIndRow, int *jacIndCol,
      blocksqp::Matrix &gradLagrange, int flag) const;
    /// Overloaded function for convenience, uses current variables of SQPiterate vars
    void calcLagrangeGradient(BlocksqpMemory* m, blocksqp::Matrix &gradLagrange, int flag) const;
    /// Update optimization tolerance (similar to SNOPT) in current iterate
    bool calcOptTol(BlocksqpMemory* m) const;

    /*
     * Solve QP subproblem
     */
    // Update the bounds on the current step, i.e. the QP variables
    void updateStepBounds(BlocksqpMemory* m, bool soc) const;
    // Solve a QP with QPOPT or qpOASES to obtain a step deltaXi and estimates
    // for the Lagrange multipliers
    int solveQP(BlocksqpMemory* m, blocksqp::Matrix &deltaXi, blocksqp::Matrix &lambdaQP,
      bool matricesChanged = true) const;
    // Compute the next Hessian in the inner loop of increasingly convexified
    // QPs and store it in vars->hess2
    void computeNextHessian(BlocksqpMemory* m, int idx, int maxQP) const;

    /*
     * Globalization Strategy
     */
    /// No globalization strategy
    int fullstep(BlocksqpMemory* m) const;
    /// Set new primal dual iterate
    void acceptStep(BlocksqpMemory* m, const blocksqp::Matrix &deltaXi,
      const blocksqp::Matrix &lambdaQP, double alpha, int nSOCS) const;
    // Overloaded function for convenience, uses current variables of SQPiterate vars
    void acceptStep(BlocksqpMemory* m, double alpha) const;
    // Reduce stepsize if a step is rejected
    void reduceStepsize(BlocksqpMemory* m, double *alpha) const;
    // Determine steplength alpha by a filter based line search similar to IPOPT
    int filterLineSearch(BlocksqpMemory* m) const;
    // Remove all entries from filter
    void initializeFilter(BlocksqpMemory* m) const;
    // Is a pair (cNorm, obj) in the current filter?
    bool pairInFilter(BlocksqpMemory* m, double cNorm, double obj) const;
    // Augment current filter by pair (cNorm, obj)
    void augmentFilter(BlocksqpMemory* m, double cNorm, double obj) const;
    // Perform a second order correction step (solve QP)
    bool secondOrderCorrection(BlocksqpMemory* m, double cNorm, double cNormTrial,
      double dfTdeltaXi, bool swCond, int it) const;
    // Reduce stepsize if a second order correction step is rejected
    void reduceSOCStepsize(BlocksqpMemory* m, double *alphaSOC) const;
    // Start feasibility restoration heuristic
    int feasibilityRestorationHeuristic(BlocksqpMemory* m) const;
    // Start feasibility restoration phase (solve NLP)
    int feasibilityRestorationPhase(BlocksqpMemory* m) const;
    // Check if full step reduces KKT error
    int kktErrorReduction(BlocksqpMemory* m) const;

    /*
     * Hessian Approximation
     */
    // Set initial Hessian: Identity matrix
    void calcInitialHessian(BlocksqpMemory* m) const;
    // [blockwise] Set initial Hessian: Identity matrix
    void calcInitialHessian(BlocksqpMemory* m, int iBlock) const;
    // Reset Hessian to identity and remove past information on Lagrange gradient and steps
    void resetHessian(BlocksqpMemory* m) const;
    // [blockwise] Reset Hessian to identity and remove past information on
    // Lagrange gradient and steps
    void resetHessian(BlocksqpMemory* m, int iBlock) const;
    // Compute current Hessian approximation by finite differences
    int calcFiniteDiffHessian(BlocksqpMemory* m) const;
    // Compute full memory Hessian approximations based on update formulas
    void calcHessianUpdate(BlocksqpMemory* m, int updateType, int hessScaling) const;
    // Compute limited memory Hessian approximations based on update formulas
    void calcHessianUpdateLimitedMemory(BlocksqpMemory* m, int updateType, int hessScaling) const;
    // [blockwise] Compute new approximation for Hessian by SR1 update
    void calcSR1(BlocksqpMemory* m, const blocksqp::Matrix &gamma, const blocksqp::Matrix &delta,
      int iBlock) const;
    // [blockwise] Compute new approximation for Hessian by BFGS update with Powell modification
    void calcBFGS(BlocksqpMemory* m, const blocksqp::Matrix &gamma, const blocksqp::Matrix &delta,
      int iBlock) const;
    // Set pointer to correct step and Lagrange gradient difference in a limited memory context
    void updateDeltaGamma(BlocksqpMemory* m) const;

    /*
     * Scaling of Hessian Approximation
     */
    // [blockwise] Update scalars for COL sizing of Hessian approximation
    void updateScalars(BlocksqpMemory* m, const blocksqp::Matrix &gamma,
      const blocksqp::Matrix &delta, int iBlock) const;
    // [blockwise] Size Hessian using SP, OL, or mean sizing factor
    void sizeInitialHessian(BlocksqpMemory* m, const blocksqp::Matrix &gamma,
      const blocksqp::Matrix &delta, int iBlock, int option) const;
    // [blockwise] Size Hessian using the COL scaling factor
    void sizeHessianCOL(BlocksqpMemory* m, const blocksqp::Matrix &gamma,
      const blocksqp::Matrix &delta, int iBlock) const;
  };

} // namespace casadi

/// \endcond
#endif // CASADI_BLOCKSQP_INTERFACE_HPP
