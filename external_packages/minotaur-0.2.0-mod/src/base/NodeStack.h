// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2008 - 2014 The MINOTAUR Team.
// 

/**
 * \file NodeStack.h
 * \brief Define the class NodeStack for storing active nodes of the
 * branch-and-bound tree on a stack data structure.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */


#ifndef MINOTAURNODESTACK_H
#define MINOTAURNODESTACK_H

#include "ActiveNodeStore.h"

namespace Minotaur {

  typedef std::deque<NodePtr> NodePtrStack;
  typedef NodePtrStack::iterator NodeStackIter;
  typedef NodePtrStack::const_iterator NodeStackConstIter;

  /** 
   * When the active nodes of the branch-and-bound tree are explored in a
   * last-in-first-out (LIFO) order, we store them in a stack. This is
   * essentially a depth first search.
   */
  class NodeStack : public ActiveNodeStore {

    public:
      /// Constructor
      NodeStack();

      /// Destroy
      virtual ~NodeStack();

      /** 
       * \brief Return true if there are no active nodes in the heap,
       * otherwise return false.
       */
      virtual bool isEmpty() const;

      /**
       * \brief Find the minimum lower bound of all the active nodes in the
       * stack.  This function is expensive and must be avoided for large
       * trees.
       */
      virtual double getBestLB() const;

      /// The maximum depth is the depth of the topmost node in the stack.
      virtual UInt getDeepestLevel() const;

      /// Remove the best node from the heap.
      virtual void pop();

      /// Write in order the node ID and the depth of each active node.
      virtual void write(std::ostream &out) const;

      /// Add a node to the set of active nodes.
      virtual void push(NodePtr n);

      /// Get access to the best node in this heap.
      virtual NodePtr top() const { return (nodes_.front()); }

      /// Get the number of active nodes in the heap.
      virtual UInt getSize() const { return (nodes_.size()); }

      /// Get iterator to the first node in the heap.
      NodeStackIter nodesBegin();

      /// Get iterator to the last node in the heap.
      NodeStackIter nodesEnd();

    private:
      /// stack of active nodes.
      NodePtrStack nodes_;
  };
  typedef boost::shared_ptr<NodeStack> NodeStackPtr;
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
