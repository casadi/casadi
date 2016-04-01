//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file NlPresHandler.h
 * \brief Declare the NlPresHandler class for try some NLP presolve ideas.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURNLPRESHANDLER_H
#define MINOTAURNLPRESHANDLER_H

#include "Handler.h"

namespace Minotaur {

class CGraph;
class CNode;
class PreAuxVars;
typedef boost::shared_ptr<CGraph> CGraphPtr;
typedef boost::shared_ptr<PreAuxVars> PreAuxVarsPtr;


/// Store statistics of presolving.
struct NlPresStats 
{
  int iters;   /// Number of iterations (main cycle).
  double time; /// Total time used in initial presolve.
  int varDel;  /// Number of variables marked for deletion.
  int conDel;  /// Number of constraints marked for deletion.
  int pRefs;   /// Number of perspective reformulations
  int vBnd;    /// Number of times variable-bounds were tightened.
  int cBnd;    /// Number of times constraint-bounds were tightened.
  int cImp;    /// Number of times coefficient in a constraint was improved.
  int nMods;   /// Number of changes in nodes.
  int qCone;   /// Number of times a quadratic constraint changed to a
               /// conic constraint.
};


/// Options for presolve.
struct NlPresOpts {
  bool doPresolve; /// True if presolve is enabled, false otherwise.
  bool showStats;  /// True if stats are displayed, false otherwise.
  int  maxIters;   /// Maximum number of iterations.
  bool coeffImp;   /// If True, do coefficient improvement.
}; 


/**
 * A NlPresHandler presolves nonlinear constraints. Experimental.
 */
class NlPresHandler : public Handler {
public:

  /// Default constructor.
  NlPresHandler();

  /// Constructor.
  NlPresHandler(EnvPtr env, ProblemPtr p);

  /// Destroy
  ~NlPresHandler();

  // Does nothing.
  void relaxInitFull(RelaxationPtr , bool *) {};

  // Does nothing.
  void relaxInitInc(RelaxationPtr , bool *) {};

  // Does nothing.
  void relaxNodeFull(NodePtr , RelaxationPtr , bool *) {};

  // Does nothing.
  void relaxNodeInc(NodePtr , RelaxationPtr , bool *) {};

  /** 
   * We assume that nonlinear presolve handler is never called for checking
   * feasibility. Always return true.
   */
  bool isFeasible(ConstSolutionPtr, RelaxationPtr, bool &, double &)
  {return true;}

  /**
   * Generate valid cuts using linear constraints.
   */
  void separate(ConstSolutionPtr , NodePtr , RelaxationPtr , CutManager *,
                SolutionPoolPtr , bool *, SeparationStatus *) {};

  /// Does nothing.
  void getBranchingCandidates(RelaxationPtr, const DoubleVector &,
                              ModVector &, BrVarCandSet &, BrCandVector &,
                              bool &) {};

  /// Does nothing.
  ModificationPtr getBrMod(BrCandPtr, DoubleVector &, 
                           RelaxationPtr, BranchDirection)
  {return ModificationPtr();}; // NULL

  /// Does nothing.
  Branches getBranches(BrCandPtr, DoubleVector &, RelaxationPtr, 
                               SolutionPoolPtr)
  {return Branches();}; // NULL

  // presolve.
  SolveStatus presolve(PreModQ *pre_mods, bool *changed);

  // Implement Handler::presolveNode().
  bool presolveNode(RelaxationPtr p, NodePtr node, SolutionPoolPtr s_pool,
                    ModVector &p_mods, ModVector &r_mods);

  // Write name
  std::string getName() const;

  /**
   * \brief Write statistics about presolve. 
   * \param [in] out The output stream to which statistics are printed.
   */
  void writePreStats(std::ostream &out) const;

  // base class method.
  void writeStats(std::ostream &out) const;
private:
  /// Should we try perspective reformulation?
  bool doPersp_;

  /// Should we try perspective reformulation?
  bool doQuadCone_;

  /// Environment.
  EnvPtr env_;

  /// Tolerance for checking feasibility etc.
  double eTol_;

  /// Log manager
  LoggerPtr logger_;
 
  /// Problem that will be presolved.
  ProblemPtr p_;

  NlPresStats stats_;

  /// Tolerance for checking zero.
  double zTol_;

  /// Who am I?
  static const std::string me_;

  void bin2Lin_(ProblemPtr p, PreModQ *mods, bool *changed);
  void bin2LinF_(ProblemPtr p, LinearFunctionPtr lf,
                 UInt nz, const UInt *irow, const UInt *jcol,
                 const double *values, PreAuxVarsPtr mods);

  bool canBin2Lin_(ProblemPtr p, UInt nz, const UInt *irow,
                   const UInt *jcol, const double *values);
  void  chkRed_(bool *changed);
  void  coeffImpr_(bool *changed);
  void  computeImpBounds_(ConstraintPtr c, VariablePtr z, 
                          double zval, double *lb, double *ub);
  void perspMod_(ConstraintPtr c, VariablePtr z);
  void perspRef_(ProblemPtr p, PreModQ *mods, bool *changed);
  void quadConeRef_(ProblemPtr p, PreModQ *mods, bool *changed);
  SolveStatus varBndsFromCons_(bool *changed);
};
typedef boost::shared_ptr<NlPresHandler> NlPresHandlerPtr;
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
