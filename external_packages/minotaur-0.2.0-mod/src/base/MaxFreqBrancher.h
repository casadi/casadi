// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file MaxFreqBrancher.h
 * \brief Declare methods and data structures for maximum-frequency branching.
 * \author Suresh B, IIT Bombay
 */


#ifndef MINOTAURMAXFREQBRANCHER_H
#define MINOTAURMAXFREQBRANCHER_H

#include "Brancher.h"

namespace Minotaur {

  class Timer;

  struct MaxFreqBrStats {
    UInt calls;      /// Number of times called to find a branching candidate.
    double time;     /// Total time spent in branching.
  };

  /// A class to select a variable for branching using maximum-frequency branching.
  class MaxFreqBrancher : public Brancher {
    
    public:
      /// Construct using an environment pointer and initialize .
      MaxFreqBrancher(EnvPtr env, HandlerVector & handlers);

      /// Destroy.
      ~MaxFreqBrancher();

      // Find a branching candidate. Returns NULL if x does not have any
      // thing to branch on or if no branching candidates are needed. 
      Branches findBranches(RelaxationPtr rel, NodePtr node, 
                            ConstSolutionPtr sol, SolutionPoolPtr s_pool, 
                            BrancherStatus & br_status, ModVector &mods);

      // Update pseudo-cost/other information after branching.
      virtual void updateAfterLP(NodePtr, ConstSolutionPtr);

      // write statistics.
      void writeStats(std::ostream &out) const;

      // Return name.
      std::string getName() const;

      /**
      * \brief Initialize data structures.
      *
      * \param[in] rel Relaxation for which this brancher is used.
      */
      void initialize(RelaxationPtr rel);

    private:
      /// Name
      const static std::string me_; 

      /// Tolerance for checking if a variable is fixed or not.
      double zTol_;

      /// Status of problem after using this brancher
      BrancherStatus status_;

      /// The problem that is being solved at this node.
      RelaxationPtr rel_;

      DoubleVector x_;

      /**
       * \brief Handlers that will be used to find candidates and the relaxations
       * after branching.
       */
      HandlerVector handlers_;

      /// Timer
      Timer *timer_;

      /// Statistics.
      MaxFreqBrStats * stats_;

      /// Set of candidates (for simple branching on variable dichotomy).
      BrVarCandSet cands_; 

      /// True if data structures are initialized. False otherwise.
      bool init_;

      /** 
       * \brief The vector that keeps a count of how many times a candidate
       * has been fraction, for each candidate.
       */
      UIntVector fracCount_;

      /** 
       * \brief The vector that keeps a count of how many times a candidate
       * has been not fixed, for each candidate.
       */
      UIntVector unfixedCount_;

      /**
       * \brief Update fractional count vector fracCount_ based upon current
       * set of candidates.
       */
      void updateFracCount_();

      /**
       * \brief Update unfixed count vector unfixedCount_ based upon current
       * relaxation.
       */
      void updateUnfixedCount_();

      /// Fills up the set of candidates in the cands_ array.
      void findCandidates_();

      /// Returns the candidate that has maximum frequency.
      BrCandPtr findBestCandidate_();

  };
  typedef boost::shared_ptr<MaxFreqBrancher> MaxFreqBrancherPtr;
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
