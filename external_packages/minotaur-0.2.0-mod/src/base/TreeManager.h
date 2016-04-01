// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

/**
 * \file TreeManager.h
 * \brief Declare class TreeManager for managing tree for Branch-and-Bound.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#ifndef MINOTAURTREEMANAGER_H
#define MINOTAURTREEMANAGER_H

#include <iostream>
#include <fstream>

#include "Types.h"

namespace Minotaur {
  
  class ActiveNodeStore;
  class WarmStart;
  typedef boost::shared_ptr<ActiveNodeStore> ActiveNodeStorePtr;
  typedef boost::shared_ptr<WarmStart> WarmStartPtr;

  // 1=like_red, 2=blue, 4=red, 5=yellow, 6=black, 7=pink, 8=cyan, 9=green
  // 11=orange, 12=green, 13=pink, 14=light blue
  /// Colors for tree-visualization using vbc
  typedef enum {
    VbcActive  = 4,   /// Unsolved, open.
    VbcFeas    = 2,   /// incumbent.
    VbcInf     = 11,  /// infeasible.
    VbcSolved  = 9,   /// solved.
    VbcSolving = 8,   /// currently being solved.
    VbcSubInf  = 13,  /// subtree is infeasible.
    VbcSubOpt  = 6    /// suboptimal.
  } VbcColors;


  /// Base class for managing the branch-and-bound tree. 
  class TreeManager {

  public:
    /// Constructor.
    TreeManager(EnvPtr env);

    /// Destroy.
    ~TreeManager();

    /// Return true if any active nodes remain in the tree. False otherwise.
    bool anyActiveNodesLeft();

    /**
     * \brief Branch and create new nodes.
     *
     * \param[in] branches The branching constraints or bounds or disjunctions
     * that are used to create the new nodes after branching.
     * \param[in] node The node that we wish to branch upon.
     * \param[in] ws The warm starting information that should be linked to
     * in the new nodes.
     * \returns The first child node.
     */
    NodePtr branch(Branches branches, NodePtr node, WarmStartPtr ws);

    /**
     * \brief Return the number of active nodes, i.e. nodes that have been
     * created, but not processed yet.
     */
    UInt getActiveNodes() const;

    /// Return the cut off value. It is INFINITY if it is not set.
    double getCutOff();

    /**
     * \brief Return the gap between the lower and upper bound as a
     * percentage. It is calculated as
     * \f$\frac{ub-lb}{\left|ub\right|+\epsilon}\times 100\f$.
     */
    double getPerGap();

    /**
     * \brief Return the value of the highest lower bound evaluated in the
     * last update of he bound.
     *
     * Since evaluating the lower bound may be expensive or certain types of
     * tree-managers, it may be updated infrequently. Consquently, the value
     * returned here may be lower than the actual lower bound of
     * the tree. Also see updateLb().
     */
    double getLb();

    /// Return the best known upper bound.
    double getUb();

    /**
     * \brief Return the size of the tree, including both active and processed
     * nodes.
     */
    UInt getSize() const;

    /**
     * \brief Search for the best candidate that can be processed next.
     *
     * It may prune some of the nodes if their lower bound is more than the
     * upper bound.  \return the best candidate found. If no candidate is
     * found, it returns NULL. The candidate is not removed from the storage.
     * It is removed only when removeActiveNode() is called.
     */
    NodePtr getCandidate();

    /**
     * \brief Insert the root node into the tree.
     *
     * \param[in] node The root node.
     */
    void insertRoot(NodePtr node);

    /**
     * \brief Prune a given node from the tree
     *
     * \param[in] node The node that must be pruned.
     */
    void pruneNode(NodePtr node);

    /**
     * \brief Remove a given active node from storage.
     *
     * It should be called after the node has been processed.
     * \param[in] node The node to be removed from storage.
     */
    void removeActiveNode(NodePtr node);

    /**
     * \brief Set the cut off value for the objective function.
     *
     * Nodes with lower bound \f$ lb \geq value-\epsilon\f$ can be pruned.
     * \param[in] value The cut off value. It can be INFINITY.
     */
    void setCutOff(double value);

    /** 
     * \brief Set the best known objective function value.
     *
     * The function does NOT check if the given value is better than the
     * previous. So care should be taken to pass only the best known value. It
     * also updates the cutoff value that is used to prune nodes.
     * \param[in] value The best known upper bound.
     */
    void setUb(double value);

    /// Return true if the tree-manager recommends diving. False otherwise.
    bool shouldDive();

    /** 
     * \brief Recalculate and return the lower bound of the tree.
     *
     * If the active nodes are stored
     * in a heap, nothing really needs to be done. For other types of storage,
     * this operation may be expensive. The result is cached into
     * bestLowerBound_.
     * \return the updated lower bound.
     */
    double updateLb();

  private:
    /// Set of nodes that are still active (those who need to be processed).
    ActiveNodeStorePtr active_nodes_; 

    /// An active node that is not in the ActiveNodeStore. One such node may
    /// exist. When we are diving. It must be deleted in the end.
    NodePtr aNode_;

    /// \brief Best known lower bound based on the last update.
    double bestLowerBound_;

    /// \brief Best known upper bound.
    double bestUpperBound_;

    /// Delete all nodes, active and inactive.
    void clearAll();

    /// The cutoff value above which nodes are assumed infeasible.
    double cutOff_;

    /// Whether we should store tree information for vbc.
    bool doVbc_;

    /// Tolerance for pruning nodes on the basis of bounds.
    const double etol_;

    /// The acceptable gap between final lb and final ub of the instance.
    const double reqGap_;

    /// The acceptable gap percentage between final lb and ub of the instance.
    double reqRelGap_;

    /// The search order: depth first, best first or something else.
    TreeSearchOrder searchType_;

    /**
     * \brief Number of nodes that have been created so far, including the
     * ones that were deleted.
     */
    UInt size_;
    
    /// Timer is used only for vbc tree emulation.
    Timer *timer_;

    /// File name to store tree information for vbc.
    std::ofstream vbcFile_;

    /// Check if the node can be pruned because of its bound.
    bool shouldPrune_(NodePtr node);

    /**
     * \brief Insert a candidate (that is not root) into the tree.
     *
     * \param[in] node The node that is to be inserted.
     * \param[in] pop_now True if  ...
     */
    void insertCandidate_(NodePtr node, bool pop_now = false);

    /**
     * \brief Remove a node from the tree.
     *
     * \param[in] node The node that is to be removed. It may remove the
     * parents of the current node also if they are no longer required.
     */
    void removeNode_(NodePtr node);
  };

  typedef boost::shared_ptr<TreeManager> TreeManagerPtr;

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
