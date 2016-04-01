//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file Node.h
 * \brief Define base class Node.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */


#ifndef MINOTAURNODE_H
#define MINOTAURNODE_H

#include "Types.h"

namespace Minotaur {

  class Node;
  class Relaxation;
  class WarmStart;
  typedef boost::shared_ptr<Node> NodePtr;
  typedef boost::shared_ptr<const Node> ConstNodePtr;
  typedef boost::shared_ptr<Relaxation> RelaxationPtr;
  typedef boost::shared_ptr<WarmStart> WarmStartPtr;
  
   /**
    * A Node is a node in the search tree or the branch-and-bound tree.
    * Associated with a node is a parent node from which the node is derived by
    * adding branching constraints. Similarly a node may have one or more
    * children each obtained by adding some branching constraints.
    * Associated with each node is an array (or vector) of modifications that
    * were applied to the parent in order to get the node. For instance, a
    * modification could be change of one or more bounds on the variables.
  
    * At each node, one can find a solution that is feasible to the original
    * instance or show that the relaxation at the node is infeasible or obtain
    * some other information. One can, at each node, solve a relaxation,
    * tighten it, preprocess it, run heuristics. 
  
    * A node that has been evaluated and does not need any further processing
    * and does not need to be branched upon is called a pruned node. A node
    * that has not been processed yet is called an active node. A node that has
    * been processed (is ready for branching or has been pruned) is called a 
    * fathomed node. A node whose all successors have been fathomed or deleted
    * can be deleted. Since we store at each node, the modification that
    * generated it from its parent, we can not delete nodes that have active
    * successors.
  
    * A Child of the current node can be obtained by applying a branch to the 
    * node. All immediate child nodes are saved in a vector.
    */
  class Node {    

  public:
    /// Default constructor.
    Node();

    /**
     * Construct a node from a parent node by applying all the modifications
     * stored in the given branch.
     */
    Node(NodePtr parentNode, BranchPtr branch);

    /// Default destructor.
    virtual ~Node();

    /// Add a child node.
    void addChild(NodePtr childNode);

    /**
     * At each node one can make several modifications to the problem.
     * Each such modification must be stored. This includes all the
     * modifications that were used to create this node from its parent
     * (while branching).
     */
    void addPMod(ModificationPtr m) { pMods_.push_back(m); }

    /**
     * At each node one can make several modifications to the relaxation.
     * Each such modification must be stored. This includes all the
     * modifications that were used to create this node from its parent
     * (while branching).
     */
    void addRMod(ModificationPtr m) { rMods_.push_back(m); }

    /**
     * Apply the modifications including the branching that were made 
     * at this node to the problem.
     */
    void applyPMods(ProblemPtr p);

    /**
     * Apply the modifications including the branching that were made 
     * at this node to the problem.
     */
    void applyRMods(RelaxationPtr rel);

    /**
     * Apply the modifications to problem and relaxation.
     */
    void applyMods(RelaxationPtr rel, ProblemPtr p);

    /// Access the children by iterating over the vector.
    NodePtrIterator childrenBegin() { return children_.begin(); }

    /// Last child.
    NodePtrIterator childrenEnd() { return children_.end(); }

    /**
     * Return the pointer to the branch that was used to create this node
     * from its parent.
     */
    BranchPtr getBranch() const { return branch_; }

    /// Return the depth of the node in the tree.
    UInt getDepth() const { return depth_; }

    /// Return the ID of this node.
    UInt getId() const { return id_; }

    /// Return the lower bound of the relaxation obtained at this node.
    double getLb() const { return lb_; }

    /// Number of children of this node.
    UInt getNumChildren() { return children_.size(); }

    /// Return a pointer to the parent node.
    NodePtr getParent() const { return parent_; }

    /// Get the status of this node.
    NodeStatus getStatus() const { return status_; }

    /// Get the tie-breaking score.
    double getTbScore() const { return tbScore_; }

    /// Get the warm start information.
    WarmStartPtr getWarmStart() { return ws_; }

    /// \todo Dont know what this is meant for.
    void makeChildOf(boost::shared_ptr<const Node> parent);

    /**
     * Get the first modification that was applied at this node to the
     * problem (and not the relaxation). Remember
     * that the vector starts with modifications that were used to branch.
     */
    ModificationConstIterator modsBegin() const { return pMods_.begin(); }

    /// End of modifications applied at this node.
    ModificationConstIterator modsEnd() const { return pMods_.end(); }

    /// Reverse iterators.
    ModificationRConstIterator modsRBegin() const { return pMods_.rbegin(); }

    /// Reverse iterators.
    ModificationRConstIterator modsREnd() const { return pMods_.rend(); }

    /// Set the status of this node.
    void setStatus(NodeStatus status) { status_ = status; }

    /// Get the tie-breaking score.
    void setTbScore(double d) { tbScore_ = d; }

    /**
     * Remove a child node from the list of children. If the node is
     * fathomed and if this list is empty, we can delete this node.
     */
    void removeChild(NodePtrIterator childNodeIter);

    /// Disassociate all children from this node.
    void removeChildren();

    /// Remove the pointer to the parent node. Used when deleting a node.
    void removeParent();

    /// Remove warm start information associated with this node.
    void removeWarmStart() { ws_.reset(); }

    /**
     * Set the ID of this node. ID of a node is unique for the given tree.
     * The treemanager must ensure this.
     */
    void setId(UInt id);

    /// Set the depth of the node in the tree.
    void setDepth(UInt depth);

    /// Set a lower bound for the relaxation at this node.
    void setLb(double value);

    /// Set warm start information
    void setWarmStart (WarmStartPtr ws);

    /**
     * Undo the modifications including the branching that were made 
     * at this node to the problem.
     */
    void undoPMods(ProblemPtr p);

    /**
     * Undo the modifications including the branching that were made 
     * at this node to the problem.
     */
    void undoRMods(RelaxationPtr r);

    /**
     * Undo the modifications including the branching that were made 
     * at this node to the problem and the relaxation.
     */
    void undoMods(RelaxationPtr rel, ProblemPtr p);

    ///Write the node
    void write(std::ostream &o) const;

  private:   
    /**
     * The branching constraints that were used to create this node from its
     * parent.
     */
    BranchPtr branch_; 

    /// Vector of immediate child nodes.
    NodePtrVector children_;

    /// Depth in the tree. Also tells how many predecessors.
    UInt depth_;

    /// Id of this node.
    UInt id_;

    /**
     * Lower bound on the relaxation at this node (not to original
     * relaxation).
     */
    double lb_;

    /**
     * Vector of modifications that are applied to the problem at this node
     * include the ones used in branching.
     */
    std::vector<ModificationPtr> pMods_;

    /**
     * Vector of modifications that are applied to the relaxation at this node
     * include the ones used in branching.
     */
    std::vector<ModificationPtr> rMods_;

    /// The parent of this node. This is NULL if the node is a root node.
    NodePtr parent_;

    /// The status of this node.
    NodeStatus status_;   

    /// Score to break tie-breaks.
    double tbScore_;

    /// The warm start information saved for this node
    WarmStartPtr ws_;

    /// Not allowed to copy a node.
    Node(const Node &node); 

    /// Not allowed to copy a node.
    Node(NodePtr node); 
  };


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
