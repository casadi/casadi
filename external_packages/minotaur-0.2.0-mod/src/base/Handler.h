//
//     MINOTAUR -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
//

/**
 * \file Handler.h
 * \brief Define abstract base class for handlers of various kinds.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 * 
 * Defines the base class Handler. 
 */

#ifndef MINOTAURHANDLER_H
#define MINOTAURHANDLER_H

#include "Types.h"

namespace Minotaur {

  class   CutManager;
  class   Node;
  class   Relaxation;
  class   PreMod;
  class   Solution;
  class   SolutionPool;
  typedef boost::shared_ptr<Node> NodePtr;
  typedef boost::shared_ptr<Relaxation> RelaxationPtr;
  typedef boost::shared_ptr<PreMod> PreModPtr;
  typedef boost::shared_ptr<const Solution> ConstSolutionPtr;
  typedef boost::shared_ptr<SolutionPool> SolutionPoolPtr;
  typedef std::deque<PreModPtr> PreModQ;
  typedef PreModQ::iterator PreModQIter;
  typedef PreModQ::const_iterator PreModQConstIter;


  /**
   * \brief Base class for handling specific types of constraints or
   * objective.
   *
   * A Handler is used to handle a particular kind of constraint/objective. 
   * It provides methods for creating a relaxation, checking if a given point
   * is feasible, separating a given point and providing a branching
   * candidate.
   */
  class Handler {
  public:

    /// Default constructor.
    Handler() {};


    /// Destroy.
    virtual ~Handler() {};

    /**
     * \brief Add constraint to be handled by this handler.
     *
     * \param[in] newcon Constraint to be added.
     */
    virtual void addConstraint(ConstraintPtr newcon)
    { cons_.push_back(newcon); };

    /**
     * \returns The beginning of constraints handled by this handler.
     */
    virtual ConstraintVector::const_iterator consBegin() const
    {return cons_.begin();};

    /**
     * \returns The end of constraints handled by this handler.
     */
    virtual ConstraintVector::const_iterator consEnd() const
    {return cons_.end();};

    /**
     * \brief Return branches for branching.
     *
     * Get branches by branching on the given candidate. In the general scheme
     * for branching, we store only bound changes, though we are also
     * capable of storing other mods. 
     *
     * \param[in] cand Candidate on which the brancher wants to branch.
     * \param[in] x The solution of the relaxation at the current node.
     * \param[in] rel The relaxation at the current node.
     * \param[in] s_pool Best feasible solutions found so far.
     * \return a vector of branch-objects.
     */
    virtual Branches getBranches(BrCandPtr cand, DoubleVector &x,
                                 RelaxationPtr rel, SolutionPoolPtr s_pool) = 0;

    /**
     * \brief find branching candidates.
     *
     * A brancher will ask each handler by calling this function to list
     * branching candidates. 
     * \param[in] rel Relaxation being solved at current node.
     * \param[in] x Solution of the relaxation.
     * \param[out] mods Any modifications that the handler found (Unused).
     * \param[out] cands The set of candidates to which branching candidates
     * must be inserted. This set has candidates that want to branch on a single
     * variable only. Other candidates must go into gencands.
     * \param[out] gencands The vector of general branching candidates. All
     * candidates that do not want to branch on a variable dichotomy must be
     * added in this vector.
     * \param[out] is_inf true if the handler finds that the problem 
     * is infeasible and the node can be pruned.
     */
    virtual void getBranchingCandidates(RelaxationPtr rel, 
                                        const DoubleVector &x,
                                        ModVector &mods, BrVarCandSet &cands,
                                        BrCandVector &gencands,
                                        bool &is_inf) = 0;

    /**
     * \brief Get the modifcation that creates a given (up or down) branch.
     *
     * If one branch is pruned by the brancher (e.g. strong brancher), 
     * then we can apply the modifications of
     * other branch to the relaxation without branching. This routine
     * returns such a modification. This function is also called for
     * obtaining modifications for strong branching on a candidate.
     * \param[in] cand The candidate for which we want the
     * modification.
     * \param[in] x The solution of relaxation at current node. 
     * \param[in] rel The relaxation at current node. 
     * \param[in] dir The Direction for which we want the modification, Up
     * or Down?.
     * \return Modification that can be applied to the relaxation before
     * re-solving it.
     */
    virtual ModificationPtr getBrMod(BrCandPtr cand, DoubleVector &x, 
                                     RelaxationPtr rel, BranchDirection dir) = 0;

    /// Return the name of the handler.
    virtual std::string getName() const = 0;

    /**
     * \brief Check if a solution is feasible.
     *
     * Check if a given solution is feasible for the constraints that are
     * handled by this handler. should_prune is true if the handler finds that
     * the problem itself is infeasible and the current node can be pruned
     * (which is different from a solution not being feasible).
     *
     * \param[in] sol The solution of the relaxation whose feasibility we
     * want to test.
     * \param[in] rel The relaxation.
     * \param[out] should_prune True if the relaxation is infeasible and we
     * can prune the node associated.
     * \param[out] inf_meas A measure of infeasibility. It may be used by
     * heuristics and other code to make the tree-search faster. Computing this
     * value is optional.
     * \return True if sol is feasible for constraints/objective asociated
     * with this handler. False if sol is not feasible.
     */
    virtual bool isFeasible(ConstSolutionPtr sol, RelaxationPtr rel,
                            bool &should_prune, double &inf_meas) = 0;

    /**
     * \brief Return true if this handler is needed for the problem.
     *
     * It is useful to know if a handler is required or not. For example, a
     * handler may be deleted if it is not needed.
     */
    virtual bool isNeeded() { return ! cons_.empty(); }

    /**
     * \brief Initial presolve.
     *
     * Do the initial presolve.  For now we will assume that presolve
     * modifies the given problem. We do not create a new problem. All
     * modifications that require post-processing for getting the solution
     * back are prepended to 'pre_mods' by the handler.
     *
     * \param[in] pre_mods A pointer to a queue of PreMod objects.
     * Modifications made by the presolver must be prepended (not appended)
     * to pre_mods. The order is important for post-solve.
     * \param[out] changed True if the presolve modified the problem.
     * \return status of presolve.
     */
    virtual SolveStatus presolve(PreModQ *pre_mods, bool *changed) = 0;

    /**
     * \brief Presolve the problem and its relaxation at a node.
     *
     * Presolve the problem and its relaxation at a given node. Bound
     * propagation and other simple modifications can be made in this function.
     * It is called after the node relaxation is setup but before it is solved.
     * Both the problem and its relaxation are presolved. Changes to the problem
     * are stored in the tree. Changes to the relaxation are optional and may or
     * may not be stored in the tree.
     * \param[in] rel Relaxation at the current node.
     * \param[in] node Current node.
     * \param[in] s_pool Pool of solutions.
     * \param[in] p_mods Unused. Modifications to the problem that must be
     * stored in this node so that they are applied to all descendant nodes as
     * well. All modifications must be appended not prepended.
     * \param[out] r_mods Modifications to the relaxation that must be stored in
     * this node so that they are applied to all descendant nodes as well.  All
     * modifications must be appended not prepended. This may be unnecessary in
     * certain algorithms.
     * \return true if Node can be pruned because infeasibility is detected.
     */
    virtual bool presolveNode(RelaxationPtr rel, NodePtr node,
                              SolutionPoolPtr s_pool, ModVector &p_mods,
                              ModVector &r_mods) = 0;

    /**
     * \brief Create root relaxation if doing full node relaxations.
     *
     * This method is used to add all the variables and constraints handled
     * by this handler, with the understanding that nodes will also be fully
     * rebuilt. The relaxation is already created, it should not be
     * freed or re-allocated.
     * \param[in,out] rel The relaxation that is being constructed.
     * \param[out] is_inf is true if the handler finds that the
     * problem is infeasible.
     */
    virtual void relaxInitFull(RelaxationPtr rel, bool *is_inf) = 0 ;

    /**
     * \brief Create root relaxation if doing incremental node relaxations.
     *
     * This method is used to add all the variables and constraints handled
     * by this handler, with the understanding that nodes will incrementally
     * relaxed. The relaxation is already created, it should not be
     * freed or re-allocated.
     * \param[in,out] rel The relaxation that is being constructed.
     * \param[out] is_inf is true if the handler finds that the
     * problem is infeasible.
     */
    virtual void relaxInitInc(RelaxationPtr rel, bool *is_inf) = 0 ;

    /**
     * \brief Create a relaxation for a node, building from scratch.
     *
     * Create a relaxation of the constraints. Either this method, or
     * relaxNodeInc should be called at each node. Here, we only make minor
     * modifications to the same relaxation.
     *
     * \param[in] node is the node for which relaxation is to be created.
     * \param[in] rel is the relaxation that is being
     * constructed. Do not allocate or re-allocate space for it. Just add
     * new variables or constraints to it.
     * \param[out] is_inf is true if the node can be pruned.
     */
    virtual void relaxNodeFull(NodePtr node, RelaxationPtr rel, bool *is_inf) = 0;

    /**
     * \brief Create an incremental relaxation for a node.
     *
     * Create a relaxation of the constraints.  Either this method, or
     * nodeRelaxFull relax should be called at root node. 
     * Usually we only make minor modifications to the same relaxation.
     * \param[in] node is the node for which relaxation is to be created.
     * \param[in] rel is the relaxation that is being constructed. Do not
     * allocate or re-allocate space for it. Just add new variables or
     * constraints to it.
     * \param[out] is_inf is true if the node can be pruned.
     */
    virtual void relaxNodeInc(NodePtr node, RelaxationPtr rel, bool *is_inf) = 0;

    /**
     * \brief add cuts to separate a given point.
     *
     * Add cuts to the relaxation to cutoff a solution. We assume that all
     * cuts are globally valid.
     *
     * \param[in] sol The solution that needs to be cut off
     * \param[in] node The node that we are currently solving.
     * \param[in] rel The relaxation at this node.
     * \param[in] cutman The CutManager where cuts should be sent.
     * \param[in] s_pool The SolutionPool containing solutions found so far.
     * \param[out] sol_found True if a new solution has been found while
     * separating
     * \param[out] status SeparationStatus returned by this routine.
     */
    virtual void separate(ConstSolutionPtr sol, NodePtr node, 
                          RelaxationPtr rel, CutManager *cutman,
                          SolutionPoolPtr s_pool, bool *sol_found,
                          SeparationStatus *status) = 0;

    /**
     * \brief Tell the handler whether the problem will be modified or the
     * relaxation will be modified or both. These modifications will be saved in
     * the tree as well.
     *
     * \param[in] mod_prob If true, modify the problem in branching and
     * presolving
     * \param[in] mod_rel If true, modify the relaxation in branching and
     * presolving.
     */
    virtual void setModFlags(bool mod_prob, bool mod_rel)
    {modProb_ = mod_prob; modRel_ = mod_rel;};

    /// Write statistics to ostream out.
    virtual void writeStats(std::ostream &) const {};

  protected:
    ConstraintVector cons_; 

    /// If true, modify the original (or transformed) problem.
    bool modProb_;

    /// If true, modify the relaxation of original (or transformed) problem.
    bool modRel_;


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
