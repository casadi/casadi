// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file ReliabilityBrancher.h
 * \brief Declare methods and data structures for reliability branching.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */


#ifndef MINOTAURRELIABILITYBRANCHER_H
#define MINOTAURRELIABILITYBRANCHER_H

#include "Brancher.h"

namespace Minotaur {

class Engine;
class Timer;
typedef boost::shared_ptr<Engine> EnginePtr;

struct RelBrStats {
  UInt bndChange;  /// Number of times variable bounds were changed.
  UInt calls;      /// Number of times called to find a branching candidate.
  UInt engProbs;   /// Number of times called to find a branching candidate.
  UInt iters;      /// Number of iterations in strong-branching.
  UInt strBrCalls; /// Number of times strong branching on a variable.
  double strTime;  /// Total time spent in strong-branching.
};


/// A class to select a variable for branching using reliability branching.
class ReliabilityBrancher : public Brancher {

public:
  /**
   * \brief Construct using an environment pointer and handlers.
   * \param [in] env Environment from which options etc. are obtained.
   * \param [in] handlers Array of handlers which are used to obtain branching
   * candidates and modifications. This array can not be changed while
   * branch-and-bound is running.
   */
  ReliabilityBrancher(EnvPtr env, HandlerVector & handlers);

  /// Destroy.
  ~ReliabilityBrancher();

  // base class function.
  Branches findBranches(RelaxationPtr rel, NodePtr node, 
                        ConstSolutionPtr sol, SolutionPoolPtr s_pool, 
                        BrancherStatus & br_status, ModVector &mods);

  /// Return value of trustCutoff parameter.
  bool getTrustCutoff();

  /// Get iteration limit of engine.
  UInt getIterLim();

  // base class function.
  std::string getName() const;

  /// Return the threshhold value.
  UInt getThresh() const;

  /**
   * \brief Initialize data structures.
   *
   * \param[in] rel Relaxation for which this brancher is used.
   */
  void initialize(RelaxationPtr rel);

  /// Set value of trustCutoff parameter.
  void setTrustCutoff(bool val);

  /**
   * \brief Set engine.
   *
   * \param[in] engine Pointer to the engine used in strong branching.
   */
  void setEngine(EnginePtr engine);

  /**
   * \brief Set iteration limit of engine.
   *
   * \param [in] k The new iteration limit.
   */
  void setIterLim(UInt k);

  /**
   * \brief Set the depth at which we stop strong branching.
   *
   * \param[in] k The new depth value.
   */
  void setMaxDepth(UInt k);

  /**
   * \brief Don't do strong branching on a cand if we did it 'k' nodes
   * or less ago.
   * \param[in] k The new minimum value.
   */
  void setMinNodeDist(UInt k);

  /**
   * \brief Set reliability threshhold
   * 
   * The pseudo cost is considered reliable if it has been updated k
   * times.
   * \param[in] k The new threshhold value.
   */
  void setThresh(UInt k);

  // base class function.
  void updateAfterLP(NodePtr node, ConstSolutionPtr sol);

  /// Write statistics.
  void writeStats(std::ostream &out) const;

private:

  /**
   * \brief Find the variable that was selected for branching.
   * 
   * This method can also find modifications based upon the
   * results of strong branching.
   * \param[in] objval Optimal objective value of the current relaxation.
   * \param[in] cutoff The cutoff value for objective function (an upper
   * bound).
   * \param[in] node The node at which we are branching.
   */
  BrCandPtr findBestCandidate_(const double objval, double cutoff, 
                               NodePtr node);

  /**
   * \brief Find and sort candidates for branching.
   *
   * Fills up the set of candidates in the cands_ array. 
   * The candidates have fractional values and
   * are sorted by the min(timesUp_[i], timesDown_[i]). The variables after
   * last_strong in the cands_ vector do not need any further strong 
   * branching.  
   */
  void findCandidates_();

  /**
   * \brief Find the score of a candidate based on its pseudo costs.
   *
   * \param[in] cand The branching candidate for which score is needed.
   * \param[in] ch_down The down score
   * \param[in] ch_up The up score.
   * \param[out] score The total score returned by this function.
   */
  void getPCScore_(BrCandPtr cand, double *ch_down, double *ch_up, 
                   double *score);

  /**
   * \brief Calculate score from the up score and down score.
   *
   * \param[in] up_score Up score.
   * \param[in] down_score Down score.
   */
  double getScore_(const double & up_score, const double & down_score);

  /**
   * \brief Check if branch can be pruned on the basis of engine status and
   * objective value.
   *
   * Check status and tell if we can prune the branch. is_rel is false if
   * the engine reports failure in convergence.
   * \param[in] chcutoff The minimum change in objective that will lead to
   * cutoff. A change greater than this value will exceed the
   * objective-cutoff.
   * \param[in] change The observed change in the objective function value.
   * \param[in] status Solve status of the engine.
   * \param[out] is_rel True if the engine status is reliable enough to use
   * the solution value as a reliable bound.
   */
  bool shouldPrune_(const double &chcutoff, const double &change, 
                    const EngineStatus & status, bool *is_rel);

  /** 
   * \brief Do strong branching on candidate.
   * \param[in] cand Candidate for strong branching.
   * \param[out] obj_up objective value estimate in up branch.
   * \param[out] obj_down objective value estimate in down branch.
   * \param[out] status_up engine status in up branch.
   * \param[out] status_down engine status in down branch.
   */
  void strongBranch_(BrCandPtr cand, double & obj_up, double & obj_down, 
                     EngineStatus & status_up, EngineStatus & status_down);

  /**
   * \brief Update Pseudocost based on the new costs.
   *
   * \param[in] i Index of the candidate.
   * \param[in] new_cost The new cost estimate.
   * \param[in] cost The vector of costs of all candidates. cost[i] is to
   * be updated
   * \param[in] count The vector that keeps a cound of how many times cost
   * has been updated for each candidate.
   */
  void updatePCost_(const int &i, const double &new_cost, 
                    DoubleVector &cost, UIntVector &count);

  /**
   * \brief Analyze the strong-branching results.
   *
   * \param[in] cand Candidate for which we performed strong branching.
   * \param[in] chcutoff The minimum change in objective function value that
   * will result in cutoff.
   * \param[in,out] change_up Change observed in the objective function value
   * in the up branch. It is modified if the engine status is error or
   * unknown.
   * \param[in,out] change_down Change observed in the objective function
   * value in the down branch. It is modified if the engine status is error or
   * unknown.
   * \param[in] status_up The engine status in up branch. 
   * \param[in] status_down The engine status in up branch.
   */
  void useStrongBranchInfo_(BrCandPtr cand, const double & chcutoff,
                            double & change_up, double & change_down, 
                            const EngineStatus & status_up,
                            const EngineStatus & status_down);

  /** 
   * \brief Display score details of the candidate.
   *
   * \param[in] cand Candidate
   * \param[in] score Score of candidate.
   * \param[in] change_up Up change.
   * \param[in] change_down Down change.
   */
  void writeScore_(BrCandPtr cand, double score, double change_up, 
                   double change_down);

  /**
   * \brief Display scores of unreliable candidates and pseudo-costs of reliable
   * ones.
   *
   * \param[in] out Outstream where scores are displayed.
   */
  void writeScores_(std::ostream &out);

  /// The engine used for strong branching.
  EnginePtr engine_;

  /// Tolerance for avoiding division by zero.
  const double eTol_;

  /**
   * \brief Handlers that will be used to find candidates and the relaxations
   * after branching.
   */
  HandlerVector handlers_;

  /// True if data structures initialized. False otherwise.
  bool init_;

  /// When did we last strong-branch on a candidate.
  UIntVector lastStrBranched_;

  /// If the depth of a node is greater than maxDepth_, then don't do any
  /// strong brancing.
  UInt maxDepth_;

  /**
   * \brief Maximum number of iterations to be performed in each call to the
   * engine.
   */
  UInt maxIterations_;

  /**
   * \brief Do not strong-branch on more than these many candidates
   */
  UInt maxStrongCands_;

  /// Name.
  const static std::string me_; 

  /**
   * \brief Don't do strong branching on a cand if we did it 'k' nodes
   * or less ago.
   */
  UInt minNodeDist_;

  /// Modifications that can be applied to the problem.
  ModVector mods_;

  /// Vector of pseudocosts for rounding down.
  DoubleVector pseudoDown_;

  /// Vector of pseudocosts for rounding up.
  DoubleVector pseudoUp_;

  /// The problem that is being solved at this node.
  RelaxationPtr rel_;

  /// A vector of candidates that have reliable pseudocosts.
  std::vector<BrCandPtr> relCands_;

  /// Statistics.
  RelBrStats * stats_;

  /// Status of problem after using this brancher.
  BrancherStatus status_;

  /// Timer to track time spent in this class.
  Timer *timer_;

  /**
   * \brief Number of times we have branched down and noted the effect on 
   * objective improvement.
   */
  UIntVector timesDown_;

  /**
   * \brief Number of times we have branched up and noted the effect on objective
   * improvement.
   */
  UIntVector timesUp_;

  /// How many times before we assume that the pseudo costs are reliable.
  UInt thresh_;

  /// Should we use the objective value to deduce a modification.
  bool trustCutoff_;

  /**
   * \brief A vector of candidates that will need strong branching. These
   * candidates will be arranged in some order of preferance.
   */
  std::vector<BrCandPtr> unrelCands_;

  /// The values of variables in the solution of the current relaxation.
  DoubleVector x_;

};
typedef boost::shared_ptr<ReliabilityBrancher> ReliabilityBrancherPtr;
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
