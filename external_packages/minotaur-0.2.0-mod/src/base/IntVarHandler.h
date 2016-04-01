//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file IntVarHandler.h
 * \brief Declare the IntVarHandler class for handling integer constrained
 * variables. It checks integrality and provides branching candidates. Does
 * not do any presolving and cut-generation.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURINTVARHANDLER_H
#define MINOTAURINTVARHANDLER_H

#include "Handler.h"

namespace Minotaur {


/**
 * IntVarHandler class considers integer variables of a problem. It only
 * checks integrality of the variables. It does not implement separate()
 * function.  Implements functions for isFeasible and branching. 
 */
class IntVarHandler : public Handler {
public:

  /// Default constructor.
  IntVarHandler();

  /// Constructor.
  IntVarHandler(EnvPtr env, ProblemPtr problem);

  /// Destroy.
  ~IntVarHandler();

  // Implement Handler::getBranches().
  Branches getBranches(BrCandPtr cand, DoubleVector & x,
                       RelaxationPtr rel, SolutionPoolPtr s_pool);
      
  // base class method
  void getBranchingCandidates(RelaxationPtr rel, 
                              const DoubleVector &x, ModVector &mods,
                              BrVarCandSet &cands, BrCandVector &gencands,
                              bool &is_inf);

  // Implement Handler::getBrMod().
  ModificationPtr getBrMod(BrCandPtr cand, DoubleVector &x, 
                           RelaxationPtr rel, BranchDirection dir);

  /// Return the integer tolerance.
  double getTol() const;

  // Write name
  std::string getName() const;

  // Check if solution is feasible.
  bool isFeasible(ConstSolutionPtr sol, RelaxationPtr relaxation, 
                  bool & should_prune, double &inf_meas);

  bool isNeeded();

  /// Does nothing.
  void postsolveGetX(const double *, UInt, DoubleVector *) {};

  /// Presolve. Don't do anything.
  SolveStatus presolve(PreModQ *, bool *) {return Finished;};

  /// Does nothing.
  bool presolveNode(RelaxationPtr, NodePtr, SolutionPoolPtr, ModVector &,
                    ModVector &)
  {return false;};

  // Does nothing.
  void relaxInitFull(RelaxationPtr rel, bool *is_inf) ;

  // Does nothing.
  void relaxInitInc(RelaxationPtr rel, bool *is_inf);

  // Does nothing.
  void relaxNodeFull(NodePtr node, RelaxationPtr rel, bool *is_inf) ;

  // Does nothing.
  void relaxNodeInc(NodePtr node, RelaxationPtr rel, bool *is_inf);

  /**
   * We do not need separation for this handler. So we just return without
   * doing anything.
   */
  void separate(ConstSolutionPtr, NodePtr , RelaxationPtr, CutManager *,
                SolutionPoolPtr, bool *, SeparationStatus *) {};

  /// Set the integer tolerance.
  void setTol(double tol);

  /// Does nothing.
  void writePreStats(std::ostream &) const {};

private:
  /// Environment.
  EnvPtr env_;

  /// True if we are doing guided dive, false otherwise.
  bool gDive_;

  /**
   * Tolerance for checking integrality.
   * If |round(x) - x| < intTol_, then it is considered to be integer
   * valued.
   */
  double intTol_;

  /// Log
  LoggerPtr logger_;

  /// For log:
  static const std::string me_;

  /// The problem for which the handler was created.
  ProblemPtr problem_;
};
typedef boost::shared_ptr<IntVarHandler> IntVarHandlerPtr;
typedef boost::shared_ptr<const IntVarHandler> ConstIntVarHandlerPtr;
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
