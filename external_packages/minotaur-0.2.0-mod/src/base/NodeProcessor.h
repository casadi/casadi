// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file NodeProcessor.h
 * \brief Define the NodeProcessor class for processing nodes in the
 * branch-and-bound algorithm.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */


#ifndef MINOTAURNODEPROCESSOR_H
#define MINOTAURNODEPROCESSOR_H

#include "Types.h"

namespace Minotaur {

  class Brancher;
  class Relaxation;
  class SolutionPool;
  class WarmStart;
  class CutManager;
  typedef boost::shared_ptr <Brancher> BrancherPtr;
  typedef boost::shared_ptr <Relaxation> RelaxationPtr;
  typedef boost::shared_ptr <SolutionPool> SolutionPoolPtr;
  typedef boost::shared_ptr <WarmStart> WarmStartPtr;
  typedef boost::shared_ptr<CutManager> CutManagerPtr;
  /**
   * NodeProcessor class is the place that does the most important tasks in
   * branch-and-bound. It processes a given node. By processing a node, we
   * mean:
   * -# presolving
   * -# solving the relaxation
   * -# checking if the relaxation or the problem itself  is infeasible
   * -# adding cuts, tightening bounds
   * -# branching decision.
   * 
   * The base class can be overridden with derived class for specific classes
   * of problems. In the default class, we do no presolving, solve the
   * relaxation by calling a continuous solver (lp_engine for LPs and
   * nlp_engine for NLPs). We do not add any cuts or tighten bounds. We branch
   * lexicographically.
   */

  class NodeProcessor {
    public:
      /// Default constructor.
      NodeProcessor() { }

      /// Destroy.
      virtual ~NodeProcessor() { }

      /// Set the brancher that will be used with this node processor.
      virtual void setBrancher(BrancherPtr brancher);

      /// Process the root node.
      virtual void processRootNode(NodePtr node, RelaxationPtr rel, 
          SolutionPoolPtr s_pool);

      /**
       * Process relaxation at the given node. If a solution is found,
       * it must be added to the solution pool. 
       */
      virtual void process(NodePtr node, RelaxationPtr rel, 
          SolutionPoolPtr s_pool) = 0;

      /// Find branches that will be used to branch at this node.
      virtual  Branches getBranches() = 0;

      /**
       * True if the node processor found at least one feasible solution while
       * processing this node.
       */
      virtual bool foundNewSolution() = 0;

      /**
       * Return the warm start information that will be used to start
       * processing children.
       */
      virtual WarmStartPtr getWarmStart() = 0;

      /// Return brancher.
      virtual BrancherPtr getBrancher() { return brancher_;};

      /// Write statistics to a given output stream
      virtual void writeStats(std::ostream &) const {};

      /// Write statistics to our own logger
      virtual void writeStats() const {};

      virtual void setCutManager(CutManager *) {};
    protected:
      /// What brancher is used for this processor
      BrancherPtr brancher_;

  };

  typedef boost::shared_ptr <NodeProcessor> NodeProcessorPtr;

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
