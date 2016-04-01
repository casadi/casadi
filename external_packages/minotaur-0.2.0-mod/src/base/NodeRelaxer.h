// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2014 The MINOTAUR Team.
// 

/**
 * \file NodeRelaxer.h
 * \brief Define the NodeRelaxer class for creating relaxation for nodes in the
 *  branch-and-bound algorithm.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */


#ifndef MINOTAURNODERELAXER_H
#define MINOTAURNODERELAXER_H

#include "Types.h"

namespace Minotaur {

class   Modification;
class   Relaxation;
typedef boost::shared_ptr<Relaxation> RelaxationPtr;

/**
 * NodeRelaxer class is used to create relaxations for nodes in a tree.
 * The relaxation at a node can be created in different ways, depending upon
 * the problem and the algorithm used. For instance, in NLP-BB, we can
 * simply relax the integrality constraints. In this case, we can just keep
 * one copy of the relaxation and change bounds in different nodes. For LP
 * based branch-and-reduce, the relaxation could be very different in a
 * node, and we may want to create new relaxation in every node.
 * 
 * Relaxations are created using handlers. Therefore, the NodeRelaxer must
 * have compatible handlers. For instance, we can not have a handler return
 * cuts in NLP-BB. Also, we can not have handler branch on a new variable
 * created in a relaxation, if we are creating in each node, new relaxations
 * from scratch.
 * 
 * This is the abstract base class.
 */
class NodeRelaxer {
public:
  /// Default constructor.
  NodeRelaxer() { }

  /// Destroy.
  virtual ~NodeRelaxer() { }

  /**Create the root node.  Prune is true if then node can be pruned.
     Also returns a vector of modifications that can be applied to the 
     root node.  This is used if "strong" (LP-based) bounds tightening is
     done. 
  */
  virtual RelaxationPtr createRootRelaxation(NodePtr rootNode, 
                                             bool &prune) = 0; 

  /**
   * Set the brancher that will be used with this node processor. dived
   * is true if this node is processed right after its parent. prune is
   * true if the node was found to be infeasible after creating the
   * relaxation.
   */
  virtual RelaxationPtr createNodeRelaxation(NodePtr node, bool dived, 
                                             bool &prune) = 0;

  /**
   * After processing the node, some node relaxers may like to make
   * changes. This function is the place to do it. diving is true if the
   * next node to be processed is a child node of the current node.
   */
  virtual void reset(NodePtr node, bool diving) = 0;

  /**
   * Return a pointer to the last relaxation that was created by this
   * relaxer.
   */
  virtual RelaxationPtr getRelaxation() = 0;
};

typedef boost::shared_ptr <NodeRelaxer> NodeRelaxerPtr;
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
