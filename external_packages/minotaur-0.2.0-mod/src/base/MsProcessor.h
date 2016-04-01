// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file MsProcessor.h
 * \brief Define multi-start node-processor for branch-and-bound
 * \author Ashutosh Mahajan, IIT Bombay
 */

#ifndef MINOTAURMSPROCESSOR_H	//renamed by pp
#define MINOTAURMSPROCESSOR_H	//renamed by pp

#include "NodeProcessor.h"

namespace Minotaur {

  class Engine;
  class Problem;
  class Solution;
  typedef boost::shared_ptr<Engine> EnginePtr;
  typedef boost::shared_ptr<const Problem> ConstProblemPtr;
  typedef boost::shared_ptr<const Solution> ConstSolutionPtr;

//  renamed by pp
  struct MBPStats {
    UInt bra;    /// Number of times relaxation became infeasible
    UInt inf;    /// Number of times relaxation became infeasible
    UInt opt;    /// Number of times relaxation gave optimal feasible solution
    UInt prob;   /// Number of times problem ocurred in solving
    UInt proc;   /// Number of nodes processed
    UInt ub;     /// Number of nodes pruned because of bound
  };

  /**
   * \brief Simple multi-start node-processor for branch-and-bound.
   *
   * MsProcessor is a derived class of NodeProcessor. It is meant to solve
   * multiple relaxations at each node with different starting points. It performs only pruning and branching in a
   * node. Does not call any presolving, cutting, or heuristic search.
   */
  class MsProcessor : public NodeProcessor {
  public:
    /// Default constructor
    MsProcessor();

    /// Constructor with a given engine.
    MsProcessor(EnvPtr env, EnginePtr engine, HandlerVector handlers_);

    /// Destroy
    ~MsProcessor();

    // True if a new solution was found while processing this node.
    bool foundNewSolution(); 

    // Find branches that will be used to branch at this node.
    Branches getBranches();

    // Get warm-start information.
    WarmStartPtr getWarmStart();

    // Implement NodeProcessor::process().
    void process(NodePtr node, RelaxationPtr rel, 
                 SolutionPoolPtr s_pool);

    // write statistics. Base class method.
    void writeStats(std::ostream &out) const; 

    // write statistics. Base class method.
    void writeStats() const; 

  protected:
    /// Branches found by this processor for this node
    Branches branches_;

    /**
     * If true, we continue to search, if engine reports error. If false,
     * we assume that the relaxation is infeasible when engine returns error.
     */
    bool contOnErr_;

    /// If lb is greater than cutOff_, we can prune this node.
    double cutOff_;

    /// Engine used to process the relaxation
    EnginePtr engine_;

    /// Status of the engine
    EngineStatus engineStatus_;

    /// All the handlers that are used for this processor
    HandlerVector handlers_;

    /// Log
    LoggerPtr logger_;

    /// For logging.
    static const std::string me_;

    /// How many new solutions were found by the processor.
    UInt numSolutions_;

    /// Relaxation that is processed by this processor.
    RelaxationPtr relaxation_;

    /// Statistics
    MBPStats stats_;

    /// Warm-start information for start processing the children
    WarmStartPtr ws_;

    /**
     * Check if the solution is feasible to the original problem. 
     * In case it is feasible, we can store the solution and update the
     * upper bound. Additionally, if the solution is optimal for the
     * current node, then the node can be pruned.
     */
    bool isFeasible_(NodePtr node, ConstSolutionPtr sol, 
                             SolutionPoolPtr s_pool, bool &should_prune);

    /// Solve the relaxation.
    void solveRelaxation_();

    /**
     * Check if a node can be pruned either because the relaxation is
     * infeasible or because the cost is too high.
     */
    bool shouldPrune_(NodePtr node, double solval, 
                              SolutionPoolPtr s_pool);

  };

  typedef boost::shared_ptr <MsProcessor> MsProcessorPtr;

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
