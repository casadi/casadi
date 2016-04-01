//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file KnapCovHandler.h
 * \brief Declare the KnapCovHandler class for handling knapsack cover 
 * constraints. It generates the cuts whenever they are needed. 
 * \author Serdar Yildiz, Argonne National Laboratory
 */

#ifndef MINOTAURPERSPCUTHANDLER_H
#define MINOTAURPERSPCUTHANDLER_H

#include "Handler.h"
#include "PerspCutGenerator.h"

namespace Minotaur {

// Pointers for handler.
class PerspCutHandler;
typedef boost::shared_ptr<PerspCutHandler> PerspCutHandlerPtr;
typedef boost::shared_ptr<const PerspCutHandler> ConstPerspCutHandlerPtr;

class Logger;
typedef boost::shared_ptr<Logger> LoggerPtr;

struct PCStats
{
  UInt perspcons;
  UInt cuts;
  double time;
};

class PerspCutHandler : public Handler {
public:

  /// Default constructor.
  PerspCutHandler();
  
  /// Constructor.
  PerspCutHandler(EnvPtr env, ProblemPtr problem);

  /// Destroy.
  ~PerspCutHandler();

  /// Does nothing.
  void relaxInitFull(RelaxationPtr, bool * ) {};
  
  /// Does nothing.
  void relaxInitInc(RelaxationPtr, bool * ) {};
  
  /// Does nothing.
  void relaxNodeFull(NodePtr, RelaxationPtr, bool * ) {};

  /// Does nothing.
  void relaxNodeInc(NodePtr, RelaxationPtr, bool * ) {};

  /// Check if solution is feasible.
  /// Checks all the constraints if they are satisfied by the given solution.
  bool isFeasible(ConstSolutionPtr sol, RelaxationPtr relaxation,
                  bool &should_prune, double &inf_meas);

  /**
   * We need separation for this handler to generate the knapsack cover
   * cuts.
   * A set of perspective cuts will be generated.
   */
  void separate(ConstSolutionPtr, NodePtr, RelaxationPtr, CutManager *cutman,
                SolutionPoolPtr, bool *, SeparationStatus * status);

  // Does nothing.
  virtual void getBranchingCandidates(RelaxationPtr,
                                      const DoubleVector &, ModVector &,
                                      BrVarCandSet &, BrCandVector &,
                                      bool &) {};

  /// Does nothing.
  virtual ModificationPtr getBrMod(BrCandPtr, DoubleVector &,
                                   RelaxationPtr, BranchDirection)
    {return ModificationPtr();};

  /// Does nothing.
  virtual Branches getBranches(BrCandPtr, DoubleVector &,
                               RelaxationPtr, SolutionPoolPtr)
    {return Branches();};

  /// Does nothing.
  SolveStatus presolve(PreModQ *, bool *) {return Finished;};

  /// Does nothing.
  virtual bool presolveNode(RelaxationPtr, NodePtr,
                            SolutionPoolPtr, ModVector &,
                            ModVector &) {return false;};
  
  /// Write name.
  std::string getName() const;

  /// Show statistics.
  void writeStats(std::ostream &) const;

  /// Return specific statistics.
  UInt PC_cuts() {return stats_->cuts;}
  double PC_time() {return stats_->time;}
  
private:
  /// Environment.
  EnvPtr env_;
  /// The problem for which the handler is created.
  ProblemPtr minlp_;
  /// Log.
  LoggerPtr logger_;
  /// Statistics.
  PCStats * stats_;
  /**
   * This is false if the current solution violates any perspective cuts.
   */
  bool isFeas_;
  /// Tolerance for accepting a new solution value: absolute threshold.
  const double solAbsTol_;
  /// Number of variables in MINLP.
  UInt numvars_;
  /// Tolerance for checking integrality.
  double intTol_;
  /// For log:
  static const std::string me_;

};

}

#endif // MINOTAURPERSPCUTHANDLER_H





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
