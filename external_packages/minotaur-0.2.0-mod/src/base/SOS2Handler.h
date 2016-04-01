//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file SOS2Handler.h
 * \brief Declare the SOS2Handler class for handling SOS type II constraints.
 * It checks integrality and provides branching candidates. Does
 * not do any presolving and cut-generation.
 * \author Ashutosh Mahajan, IIT Bombay
 */

#ifndef MINOTAURSOS2HANDLER_H
#define MINOTAURSOS2HANDLER_H

#include "Handler.h"

namespace Minotaur {

/**
 * SOS2Handler class considers all SOS Type 2 constraints in the problem. 
 * It does not implement separate() function.  Implements functions for
 * checking feasibility and branching. 
 */
class SOS2Handler : public Handler {
public:

  /// Default constructor.
  SOS2Handler();

  /// Constructor.
  SOS2Handler(EnvPtr env, ProblemPtr problem);

  /// Destroy.
  ~SOS2Handler();

  // Does nothing.
  void relaxInitFull(RelaxationPtr rel, bool *is_inf);

  // Does nothing.
  void relaxInitInc(RelaxationPtr rel, bool *is_inf);

  // Does nothing.
  void relaxNodeFull(NodePtr node, RelaxationPtr rel, bool *is_inf);

  // Does nothing.
  void relaxNodeInc(NodePtr node, RelaxationPtr rel, bool *is_inf);

  // Check if solution is feasible.
  bool isFeasible(ConstSolutionPtr sol, RelaxationPtr relaxation, 
                  bool &should_prune, double &inf_meas);

  /// Empty method as we do not do any separation in SOS1.
  void separate(ConstSolutionPtr, NodePtr , RelaxationPtr, CutManager *,
                SolutionPoolPtr, bool *, SeparationStatus *);

  // Base class method.
  void getBranchingCandidates(RelaxationPtr rel, const DoubleVector &x,
                              ModVector & mods, BrVarCandSet &cands,
                              BrCandVector &gencands, bool & is_inf);

  // Base class method.
  ModificationPtr getBrMod(BrCandPtr cand, DoubleVector &x, 
                           RelaxationPtr rel, BranchDirection dir);

  //  Base class method.
  Branches getBranches(BrCandPtr cand, DoubleVector & x,
                       RelaxationPtr rel, SolutionPoolPtr s_pool);
      
  bool isNeeded();

  /// Presolve. Do not do any presolving.
  SolveStatus presolve(PreModQ *, bool *) {return Finished;};

  /// Does nothing.
  void postsolveGetX(const double *, UInt, DoubleVector *);

  /// Does nothing.
  bool presolveNode(RelaxationPtr, NodePtr, SolutionPoolPtr, ModVector &,
                    ModVector &)
  {return false;};

  // Write name
  std::string getName() const;

  void getNzNumSum_(SOSPtr sos, const DoubleVector x, int *nz,
                    double *nzsum);

  /// Does nothing.
  void writePreStats(std::ostream &) const {};

  /// Return the integer tolerance.
  double getTol() const;

  /// Set the integer tolerance.
  void setTol(double tol);

private:
  /// Environment.
  EnvPtr env_;

  /// Log
  LoggerPtr logger_;

  /// For log:
  static const std::string me_;

  /// The problem for which the handler was created.
  ProblemPtr problem_;

  /**
   * \brief Tolerance for checking integrality.
   *
   * If abs(x) < zTol_, then it is considered to be zero.
   */
  double zTol_;

  void getNzAvgWt_(SOSPtr sos, const DoubleVector x,
                   const double *weights, int *nz, double *avgwt);

  void getSumN_(SOSPtr sos, const DoubleVector &x, double *sum, int *nnz,
                int *nv, int *nspos);

  bool isXFeasible_(const double *x, SOSPtr sos);

};
typedef boost::shared_ptr<SOS2Handler> SOS2HandlerPtr;
typedef boost::shared_ptr<const SOS2Handler> ConstSOS2HandlerPtr;
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
