//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file LinearHandler.h
 * \brief Declare the LinearHandler class for handling integer and continuous
 * variables.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURLINEARHANDLER_H
#define MINOTAURLINEARHANDLER_H

#include "Handler.h"

namespace Minotaur {

class LinearFunction;
typedef boost::shared_ptr<LinearFunction> LinearFunctionPtr;

/// Store statistics of presolving.
struct LinPresolveStats 
{
  int iters;   ///> Number of iterations (main cycle).
  double time; ///> Total time used in initial presolve.
  double timeN;///> Total time used in presolveNode.
  int varDel;  ///> Number of variables marked for deletion.
  int conDel;  ///> Number of constraints marked for deletion.
  int var2Bin; ///> Number of variables converted to binary.
  int var2Int; ///> Number of variables converted to integers.
  int vBnd;    ///> Number of times variable-bounds were tightened.
  int cBnd;    ///> Number of times constraint-bounds were tightened.
  int cImp;    ///> Number of times coefficient in a constraint was improved.
  int bImpl;   ///> No. of times a binary var. was changed to implied binary.
  int nMods;   ///> Number of changes made in all nodes.
};

/// Options for presolve.
struct LinPresolveOpts {
  bool doPresolve; /// True if presolve is enabled, false otherwise.

  bool showStats;  /// True if stats are displayed, false otherwise.

  int  maxIters;   /// Maximum number of iterations.

  bool purgeVars;  /// If True, purge fixed variables.

  bool purgeCons;  /// If True, purge redundant constraints.

  bool dualFix;    /// If True, do dual cost fixing.

  bool coeffImp;   /// If True, do coefficient improvement.
}; 


/**
 * An LinearHandler handles variables of a problem. It only checks bounds 
 * and integrality of the variables.
 */
class LinearHandler : public Handler {
public:

  /// Default constructor.
  LinearHandler();

  /// Constructor.
  LinearHandler(EnvPtr env, ProblemPtr problem);

  /// Destroy
  ~LinearHandler();

  // Does nothing.
  void relaxInitFull(RelaxationPtr rel, bool *is_inf) ;

  // Does nothing.
  void relaxInitInc(RelaxationPtr rel, bool *is_inf);

  // Does nothing.
  void relaxNodeFull(NodePtr node, RelaxationPtr rel, bool *is_inf) ;

  // Does nothing.
  void relaxNodeInc(NodePtr node, RelaxationPtr rel, bool *is_inf);

  /** 
   * We assume that linear constraints and bound constraints are always
   * satisfied. Always return true.
   */
  bool isFeasible(ConstSolutionPtr, RelaxationPtr, bool &, double &)
  {return true;}

  bool isNeeded() { return true; }

  /**
   * Generate valid cuts using linear constraints.
   */
  void separate(ConstSolutionPtr sol, NodePtr node, 
                RelaxationPtr rel, CutManager *cutman, SolutionPoolPtr s_pool, 
                bool *sol_found, SeparationStatus *status);

  /// Does nothing.
  virtual void getBranchingCandidates(RelaxationPtr , 
                                      const DoubleVector &, ModVector &, 
                                      BrVarCandSet &, BrCandVector &,
                                      bool &) {};

  /// Does nothing.
  virtual ModificationPtr getBrMod(BrCandPtr, DoubleVector &, 
                                   RelaxationPtr, BranchDirection)
  {return ModificationPtr();}; // NULL

  /// Does nothing.
  virtual Branches getBranches(BrCandPtr, DoubleVector &, RelaxationPtr, 
                               SolutionPoolPtr)
  {return Branches();}; // NULL

  // presolve.
  virtual SolveStatus presolve(PreModQ *pre_mods, bool *changed);

  // Implement Handler::presolveNode().
  virtual bool presolveNode(RelaxationPtr p, NodePtr node,
                            SolutionPoolPtr s_pool, ModVector &p_mods,
                            ModVector &r_mods);
      

  // Write name
  virtual std::string getName() const;

  /// Return a constant pointer to the presolve options.
  const LinPresolveOpts* getOpts() const;

  /// If true show statistics.
  void setPreOptShowStats(bool val) {pOpts_->showStats = val;}; 

  /// Maximum number of iterations.
  void setPreOptMaxIters(int val) {pOpts_->maxIters = val;}; 

  /// If True, purge fixed variables.
  void setPreOptPurgeVars(bool val) {pOpts_->purgeVars = val;}; 

  /// If True, purge redundant constraints.
  void setPreOptPurgeCons(bool val) {pOpts_->purgeCons = val;}; 

  void setPreOptDualFix(bool val) {pOpts_->dualFix = val;}; 

  void setPreOptCoeffImp(bool val) {pOpts_->coeffImp = val;}; 

  void simplePresolve(ProblemPtr p, SolutionPoolPtr spool, ModVector &t_mods,
                      SolveStatus &status);

  /// Write the presolve statistics.
  void writePreStats(std::ostream &out) const;

  // Write statistics.
  void writeStats(std::ostream &out) const;

protected:
  /// Environment.
  EnvPtr env_;

  /// The problem for which the handler was created.
  ProblemPtr problem_;

  /// Log
  LoggerPtr logger_;

  /// Tolerance for checking integrality.
  /**
   * If |round(x) - x| < intTol_, then it is considered to be integer
   * valued.
   */
  const double intTol_;

  /// Tolerance.
  const double eTol_;

  /// Infinity. Bounds beyond this number are treated as infinity.
  const double infty_;

  /// Statistics of presolve.
  LinPresolveStats *pStats_;

  /// Options for presolve.
  LinPresolveOpts *pOpts_;

  /**
   * Linear variables: variables that do not appear in nonlinear
   * functions, both in objective and constraints.
   */
  VarQueue linVars_;


  /// For log.
  static const std::string me_;

  void chkIntToBin_(VariablePtr v);

  void chkSing_(bool *changed);
  void coeffImp_(bool *changed);
  void computeImpBounds_(ConstraintPtr c, VariablePtr z, double zval,
                         double *lb, double *ub);
  void copyBndsFromRel_(RelaxationPtr rel, ModVector &p_mods);

  void delFixedVars_(bool *changed);

  void dualFix_(bool *changed);
  void dupRows_(bool *changed);

  /// check if lb <= ub for all variables and constraints.
  SolveStatus checkBounds_(ProblemPtr p);

  void findLinVars_();

  void findAllBinCons_();
  void fixToCont_();

  void getLfBnds_(LinearFunctionPtr lf, double *lo, double *up);
  void getSingLfBnds_(LinearFunctionPtr lf, double *lo, double *up);

  SolveStatus linBndTighten_(ProblemPtr p, bool apply_to_prob, 
                      ConstraintPtr c_ptr, bool *changed, ModQ *mods, UInt *nintmods);

  void purgeVars_(PreModQ *pre_mods);

  /**
   * \brief Common routine for building relaxation by copying all the linear
   * constraints and variable-bounds from a given problem.
   * 
   * \param[in] p The problem whose relaxation we want to create.
   * \param[in] rel The relaxation in which we want to add new variables
   * and constraints.
   * \param [out] is_inf True if problem p is found to be infeasible, false
   * otherwise.
   */
  void relax_(ProblemPtr p, RelaxationPtr rel, bool *is_inf);

  void substVars_(bool *changed, PreModQ *pre_mods);

  /// Round the bounds
  void tightenInts_(ProblemPtr p, bool apply_to_prob, bool *changed, 
                    ModQ *mods);

  bool treatDupRows_(ConstraintPtr c1, ConstraintPtr c2, double mult,
                     bool *changed);

  void updateLfBoundsFromLb_(ProblemPtr p, bool apply_to_prob, 
                             LinearFunctionPtr lf, double lb, double uu,
                             bool is_sing, bool *changed, ModQ *mods,
                             UInt *nintmods);

  void updateLfBoundsFromUb_(ProblemPtr p, bool apply_to_prob, 
                             LinearFunctionPtr lf, double ub, double ll,
                             bool is_sing, bool *changed, ModQ *mods,
                             UInt *nintmods);

  SolveStatus varBndsFromCons_(ProblemPtr p, bool apply_to_prob, bool *changed, 
                               ModQ *mods, UInt *nintmods);

  SolveStatus varBndsFromObj_(ProblemPtr p, double ub, bool apply_to_prob, 
                              bool *changed, ModQ *mods);
};
typedef boost::shared_ptr<LinearHandler> LinearHandlerPtr;
typedef boost::shared_ptr<const LinearHandler> ConstLinearHandlerPtr;
}
#endif

// Local Variables: 
// mode: c++ 
// eval: (c-set-style "k&r") 
// eval: (c-set-offset 'innamespace 0) 
// eval: (setq c-basic-offset 2) 
// eval: (setq fill-column 78) 
// eval: (auto-fill-mode 1) 
// eval: (setq column-number-mode 1) 
// eval: (setq indent-tabs-mode nil) 
// End:
